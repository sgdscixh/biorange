import os
import pandas as pd
import gzip
import mygene
from typing import Generator, Union, List, Dict
import tempfile
from biorange.logger import get_logger
from biorange.utils.package_fileload import get_data_file_path

# 设置日志记录
logger = get_logger(__name__)


class TCMDataProcessor:
    def __init__(self):
        self.mg = mygene.MyGeneInfo()
        self.temp_dir = tempfile.TemporaryDirectory()
        logger.info(f"Temporary directory created at {self.temp_dir.name}")

    def _read_csv(self, file_path: str, sep: str = ",") -> pd.DataFrame:
        """读取 CSV 文件并处理可能的空文件错误。"""
        try:
            df = pd.read_csv(file_path, sep=sep)
            if df.empty:
                logger.warning(f"The file {file_path} is empty.")
            return df
        except pd.errors.EmptyDataError:
            logger.error(
                f"No data found in {file_path}. Please check the file content."
            )
            return pd.DataFrame()
        except Exception as e:
            logger.error(f"Error reading {file_path}: {e}")
            return pd.DataFrame()

    def _read_large_gzipped_tsv(
        self, file_path: str, chunksize: int = 100000
    ) -> Generator[pd.DataFrame, None, None]:
        """读取大型 gzipped TSV 文件，并处理可能的空文件错误。"""
        try:
            with gzip.open(file_path, "rt") as f:
                for chunk in pd.read_csv(f, sep="\t", chunksize=chunksize):
                    yield chunk
        except pd.errors.EmptyDataError:
            logger.error(
                f"No data found in {file_path}. Please check the file content."
            )
        except Exception as e:
            logger.error(f"Error reading {file_path}: {e}")

    def _merge_files(
        self, df1: pd.DataFrame, df2: pd.DataFrame, on: str, how: str = "left"
    ) -> pd.DataFrame:
        """通用合并函数，处理空文件情况。"""
        if df1.empty or df2.empty:
            logger.error("One or more input dataframes are empty, skipping merge.")
            return pd.DataFrame()
        return pd.merge(df1, df2, on=on, how=how)

    def _save_dataframe(self, df: pd.DataFrame, filename: str) -> str:
        """保存 DataFrame 为 CSV 文件到临时目录，并返回文件路径。"""
        temp_file_path = os.path.join(self.temp_dir.name, filename)
        df.to_csv(temp_file_path, index=False)
        logger.info(f"Data saved to {temp_file_path}")
        return temp_file_path

    def _merge_with_large_chemical_data(
        self, inchikeys: Union[str, List[str]], gzipped_file: str
    ) -> pd.DataFrame:
        """将 InChIKey 数据与大型化学物质数据进行匹配合并。"""
        if isinstance(inchikeys, str):
            inchikey_list = [inchikeys]
        else:
            inchikey_list = inchikeys

        input_df = pd.DataFrame({"inchikey": inchikey_list})
        if input_df.empty:
            logger.error(f"InChIKey list is empty. Skipping this merge step.")
            return pd.DataFrame()

        results = [
            pd.merge(
                input_df,
                chunk[["flat_chemical_id", "inchikey"]],
                on="inchikey",
                how="inner",
            )
            for chunk in self._read_large_gzipped_tsv(gzipped_file)
        ]

        result_df = pd.concat(results, ignore_index=True)
        return result_df

    def _map_chemical_to_protein(
        self, chemical_df: pd.DataFrame, protein_file: str
    ) -> pd.DataFrame:
        """将化学物质 ID 映射到蛋白质，并输出中间结果。"""
        if chemical_df.empty:
            logger.error(f"Chemical DataFrame is empty. Skipping this mapping step.")
            return pd.DataFrame()

        results = [
            chemical_df.merge(
                chunk[["chemical", "protein", "combined_score"]],
                left_on="flat_chemical_id",
                right_on="chemical",
                how="left",
            )
            for chunk in self._read_large_gzipped_tsv(protein_file)
        ]

        result_df = pd.concat(results, ignore_index=True)
        return result_df

    def _convert_protein_to_gene_names(self, protein_df: pd.DataFrame) -> pd.DataFrame:
        """将 protein ID 转换为基因名，并包含 combined_score 列。"""
        if protein_df.empty:
            logger.error(f"Protein DataFrame is empty. Skipping gene name conversion.")
            return pd.DataFrame()

        protein_df["ENSP"] = protein_df["protein"].str.split(".").str[1]
        ensembl_ids = protein_df["ENSP"].tolist()

        gene_info = self.mg.querymany(
            ensembl_ids, scopes="ensembl.protein", fields="symbol", species="human"
        )
        results = [
            {"ENSP": gene["query"], "gene_name": gene.get("symbol", "N/A")}
            for gene in gene_info
        ]

        gene_df = pd.DataFrame(results)
        final_df = protein_df.merge(gene_df, on="ENSP", how="left")
        final_df["source"] = "STITCH"
        final_df = final_df[
            ["inchikey", "ENSP", "gene_name", "combined_score", "source"]
        ]

        # 删除空白行
        final_df.dropna(subset=["gene_name", "combined_score"], inplace=True)

        # 添加调试信息
        logger.info(
            f"Number of rows with missing gene_name: {final_df['gene_name'].isnull().sum()}"
        )
        logger.info(
            f"Number of rows with missing combined_score: {final_df['combined_score'].isnull().sum()}"
        )

        return final_df

    def search(
        self,
        inchikeys: Union[str, List[str]],
        chemical_file: str = get_data_file_path("TCMSP_NGM_STITCH_INCHIKEY.tsv.gz"),
        protein_file: str = get_data_file_path(
            "9606.protein_chemical.links.transfer.v5.0.tsv.gz"
        ),
    ) -> pd.DataFrame:
        """主接口：根据 InChIKey 查找对应的基因名。"""

        chemical_df = self._merge_with_large_chemical_data(inchikeys, chemical_file)
        if chemical_df.empty:
            logger.error(
                "Failed to merge InChIKey with chemical data. Aborting search."
            )
            return pd.DataFrame()

        protein_df = self._map_chemical_to_protein(chemical_df, protein_file)
        if protein_df.empty:
            logger.error("Failed to map chemicals to proteins. Aborting search.")
            return pd.DataFrame()

        final_df = self._convert_protein_to_gene_names(protein_df)
        if final_df.empty:
            logger.error("Failed to convert proteins to gene names.")
            return pd.DataFrame()

        # 过滤后的 DataFrame
        filtered_df = final_df[final_df["combined_score"] > 300]

        return filtered_df

    def get_rawdata(
        self,
        inchikeys: Union[str, List[str]],
        chemical_file: str = get_data_file_path("TCMSP_NGM_STITCH_INCHIKEY.tsv.gz"),
        protein_file: str = get_data_file_path(
            "9606.protein_chemical.links.transfer.v5.0.tsv.gz"
        ),
    ) -> pd.DataFrame:
        """获取原始数据并保存到文件。"""
        chemical_df = self._merge_with_large_chemical_data(inchikeys, chemical_file)
        if chemical_df.empty:
            logger.error(
                "Failed to merge InChIKey with chemical data. Aborting search."
            )
            return pd.DataFrame()

        protein_df = self._map_chemical_to_protein(chemical_df, protein_file)
        if protein_df.empty:
            logger.error("Failed to map chemicals to proteins. Aborting search.")
            return pd.DataFrame()

        final_df = self._convert_protein_to_gene_names(protein_df)
        if final_df.empty:
            logger.error("Failed to convert proteins to gene names.")
            return pd.DataFrame()

        # 将最终结果保存到指定路径
        final_file_path = self._save_dataframe(final_df, "stitch_target_raw.csv")
        logger.info(f"Final raw data saved to {final_file_path}")

        return final_df


stich_inchikey_target = TCMDataProcessor().search

# 使用示例
if __name__ == "__main__":
    processor = TCMDataProcessor()
    df = pd.read_csv(
        "/home/liuyan/projects/package/biorange/biorange/target_predict/input_data/admet_filtered_ingredients.csv"
    )
    df_list = df["inchikey"].tolist()
    result_df = processor.search(
        inchikeys=df_list,
    )
    result_df = result_df.drop_duplicates()  # 去除重复行
    print(result_df)
    # result_df.to_csv("results/output2/stitch_target_raw333.csv", index=False)
