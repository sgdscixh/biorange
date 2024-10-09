import json
import re
import time
import zlib
from typing import Any, Dict, Generator, List, Optional, Union
from urllib.parse import parse_qs, urlencode, urlparse
from xml.etree import ElementTree
import logging

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry


logger = logging.getLogger(__name__)

# Constants
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

# Retry strategy for requests
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


class UniportMapping:
    def __init__(self):
        self.session = session

    def check_response(self, response: requests.Response) -> None:
        try:
            response.raise_for_status()
        except requests.HTTPError as e:
            logger.error(f"HTTPError: {e.response.json()}")
            raise

    def submit_id_mapping(self, from_db: str, to_db: str, ids: List[str]) -> str:
        response = self.session.post(
            f"{API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
            timeout=60000,
        )
        self.check_response(response)
        return response.json()["jobId"]

    def get_next_link(self, headers: Dict[str, str]) -> Optional[str]:
        if "Link" in headers:
            match = re.match(r'<(.+)>; rel="next"', headers["Link"])
            return match.group(1) if match else None
        return None

    def check_id_mapping_results_ready(self, job_id: str) -> bool:
        while True:
            response = self.session.get(f"{API_URL}/idmapping/status/{job_id}")
            self.check_response(response)
            status = response.json().get("jobStatus")
            if status == "RUNNING":
                logger.info(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            elif status:
                raise Exception(status)
            else:
                return bool(
                    response.json().get("results") or response.json().get("failedIds")
                )

    def get_batch(
        self, batch_url: str, file_format: str, compressed: bool
    ) -> Generator[Union[Dict[str, Any], List[str]], None, None]:
        while batch_url:
            response = self.session.get(batch_url)
            self.check_response(response)
            yield self.decode_results(response, file_format, compressed)
            batch_url = self.get_next_link(response.headers)

    def combine_batches(
        self,
        all_results: Union[Dict[str, Any], List[str]],
        batch_results: Union[Dict[str, Any], List[str]],
        file_format: str,
    ) -> Union[Dict[str, Any], List[str]]:
        if file_format == "json":
            for key in ("results", "failedIds"):
                if key in batch_results:
                    all_results[key].extend(batch_results[key])
        elif file_format in {"tsv", "xml"}:
            all_results.extend(
                batch_results[1:] if file_format == "tsv" else batch_results
            )
        else:
            all_results += batch_results
        return all_results

    def get_id_mapping_results_link(self, job_id: str) -> str:
        response = self.session.get(f"{API_URL}/idmapping/details/{job_id}")
        self.check_response(response)
        return response.json()["redirectURL"]

    def decode_results(
        self, response: requests.Response, file_format: str, compressed: bool
    ) -> Union[Dict[str, Any], List[str], str]:
        content = (
            zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
            if compressed
            else response.content
        )
        if file_format == "json":
            return json.loads(content.decode("utf-8"))
        elif file_format == "tsv":
            return content.decode("utf-8").splitlines()
        elif file_format in {"xlsx", "xml"}:
            return [content]
        return content.decode("utf-8")

    def get_xml_namespace(self, element: ElementTree.Element) -> str:
        match = re.match(r"\{(.*)\}", element.tag)
        return match.group(1) if match else ""

    def merge_xml_results(self, xml_results: List[str]) -> str:
        merged_root = ElementTree.fromstring(xml_results[0])
        namespace = self.get_xml_namespace(merged_root[0])
        for result in xml_results[1:]:
            root = ElementTree.fromstring(result)
            for child in root.findall(f"{{{namespace}}}entry"):
                merged_root.append(child)
        ElementTree.register_namespace("", namespace)
        return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)

    def print_progress_batches(self, batch_index: int, size: int, total: int) -> None:
        n_fetched = min((batch_index + 1) * size, total)
        logger.info(f"Fetched: {n_fetched} / {total}")

    def get_id_mapping_results_search(self, url: str) -> Union[Dict[str, Any], str]:
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query.get("format", ["json"])[0]

        # Ensure 'size' is set to 500 if not present
        if "size" not in query:
            query["size"] = ["500"]
        size = int(query["size"][0])

        compressed = query.get("compressed", ["false"])[0].lower() == "true"
        parsed = parsed._replace(query=urlencode(query, doseq=True))
        url = parsed.geturl()

        response = self.session.get(url)
        self.check_response(response)
        results = self.decode_results(response, file_format, compressed)
        total = int(response.headers["x-total-results"])
        self.print_progress_batches(0, size, total)

        for i, batch in enumerate(
            self.get_batch(
                self.get_next_link(response.headers), file_format, compressed
            ),
            1,
        ):
            results = self.combine_batches(results, batch, file_format)
            self.print_progress_batches(i, size, total)

        return self.merge_xml_results(results) if file_format == "xml" else results

    def get_id_mapping_results_stream(
        self, url: str
    ) -> Union[Dict[str, Any], List[str]]:
        if "/stream/" not in url:
            url = url.replace("/results/", "/results/stream/")
        response = self.session.get(url)
        self.check_response(response)
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query.get("format", ["json"])[0]
        compressed = query.get("compressed", ["false"])[0].lower() == "true"
        return self.decode_results(response, file_format, compressed)

    # def convert_results_to_dataframe(self, results: Dict[str, Any]) -> pd.DataFrame:
    #     rows = [
    #         [
    #             result["from"],
    #             result["to"]["primaryAccession"],
    #             (
    #                 result["to"]["genes"][0]["geneName"]["value"]
    #                 if result["to"]["genes"]
    #                 else None
    #             ),
    #             result["to"]["organism"]["scientificName"],
    #         ]
    #         for result in results["results"]
    #     ]
    #     return pd.DataFrame(
    #         rows, columns=["chembal", "uniport_accession", "gene_name", "organism"]
    #     )
    def convert_results_to_dataframe(self, results: Dict[str, Any]) -> pd.DataFrame:
        rows = []
        for result in results["results"]:
            gene_name = None
            if "genes" in result["to"] and result["to"]["genes"]:
                gene_name = result["to"]["genes"][0].get("geneName", {}).get("value")

            rows.append(
                [
                    result["from"],
                    result["to"]["primaryAccession"],
                    gene_name,
                    result["to"]["organism"]["scientificName"],
                ]
            )

        return pd.DataFrame(
            rows, columns=["chembal", "uniport_accession", "gene_name", "organism"]
        )

    def get_dataframe_from_ids(self, ids: List[str] | str) -> pd.DataFrame:
        job_id = self.submit_id_mapping(from_db="ChEMBL", to_db="UniProtKB", ids=ids)
        if self.check_id_mapping_results_ready(job_id):
            link = self.get_id_mapping_results_link(job_id)
            results_dict = self.get_id_mapping_results_search(link)
            return self.convert_results_to_dataframe(results_dict)
        return pd.DataFrame()


if __name__ == "__main__":

    def fetch_and_merge_data(df):
        scraper = UniportMapping()
        res = scraper.get_dataframe_from_ids(df["ids"])

        # 合并 DataFrame
        merged_query = pd.merge(
            df,
            res,
            left_on="ids",
            right_on="chembal",
            how="left",
        )
        return merged_query

    data = {
        "ids": ["CHEMBL2902", "CHEMBL2725", "CHEMBL1234"],
        "other_info": ["info1", "info2", "info3"],
    }
    df = pd.DataFrame(data)

    # 调用函数获取合并后的 DataFrame
    merged_df = fetch_and_merge_data(df)

    # 打印合并后的结果
    print(merged_df)
