import os
import pandas as pd
import gzip
import mygene
from typing import Generator
import tempfile
from biorange.logger import get_logger

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
        self, inchikey_file: str, gzipped_file: str
    ) -> pd.DataFrame:
        """将 InChIKey 数据与大型化学物质数据进行匹配合并。"""
        input_df = self._read_csv(inchikey_file)
        if input_df.empty:
            logger.error(f"{inchikey_file} is empty. Skipping this merge step.")
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
        chemical_file_path = self._save_dataframe(
            result_df, "TCM_stitch_chemical_id.csv"
        )
        return result_df, chemical_file_path

    def _map_chemical_to_protein(
        self, chemical_file: str, protein_file: str
    ) -> pd.DataFrame:
        """将化学物质 ID 映射到蛋白质，并输出中间结果。"""
        input_df = self._read_csv(chemical_file)
        if input_df.empty:
            logger.error(f"{chemical_file} is empty. Skipping this mapping step.")
            return pd.DataFrame()

        results = [
            input_df.merge(
                chunk[["chemical", "protein", "combined_score"]],
                left_on="flat_chemical_id",
                right_on="chemical",
                how="left",
            )
            for chunk in self._read_large_gzipped_tsv(protein_file)
        ]

        result_df = pd.concat(results, ignore_index=True)
        protein_file_path = self._save_dataframe(result_df, "stitch_ensp_score.csv")
        return result_df, protein_file_path

    def _convert_protein_to_gene_names(self, protein_file: str) -> pd.DataFrame:
        """将 protein ID 转换为基因名，并包含 combined_score 列。"""
        df = self._read_csv(protein_file)
        if df.empty:
            logger.error(f"{protein_file} is empty. Skipping gene name conversion.")
            return pd.DataFrame()

        df["ENSP"] = df["protein"].str.split(".").str[1]
        ensembl_ids = df["ENSP"].tolist()

        gene_info = self.mg.querymany(
            ensembl_ids, scopes="ensembl.protein", fields="symbol", species="human"
        )
        results = [
            {"ENSP": gene["query"], "gene_name": gene.get("symbol", "N/A")}
            for gene in gene_info
        ]

        gene_df = pd.DataFrame(results)
        final_df = df.merge(gene_df, on="ENSP", how="left")
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

        final_file_path = self._save_dataframe(final_df, "stitch_target_raw.csv")

        # 生成过滤后的文件
        filtered_df = final_df[final_df["combined_score"] > 300]
        filtered_file_path = self._save_dataframe(
            filtered_df, "stitch_output_target.csv"
        )

        # 提取未查询到基因名的ENSP
        ensp_without_gene_name_df = df[
            df["ENSP"].isin(gene_df[gene_df["gene_name"] == "N/A"]["ENSP"])
        ]
        ensp_without_gene_name_file_path = self._save_dataframe(
            ensp_without_gene_name_df[["inchikey", "ENSP", "combined_score"]],
            "ENSP_without_gene_name.csv",
        )

        return final_df, filtered_df, ensp_without_gene_name_df

    def search(
        self,
        inchikey_file: str,
        chemical_file: str,
        protein_file: str,
        output_raw_file: str = "stitch_target_raw.csv",
        output_filtered_file: str = "stitch_output_target.csv",
        output_ensp_without_gene_name_file: str = "ENSP_without_gene_name.csv",
    ) -> pd.DataFrame:
        """主接口：根据 InChIKey 查找对应的基因名。"""
        chemical_df, chemical_file_path = self._merge_with_large_chemical_data(
            inchikey_file, chemical_file
        )
        if chemical_df.empty:
            logger.error(
                "Failed to merge InChIKey with chemical data. Aborting search."
            )
            return pd.DataFrame()

        protein_df, protein_file_path = self._map_chemical_to_protein(
            chemical_file_path, protein_file
        )
        if protein_df.empty:
            logger.error("Failed to map chemicals to proteins. Aborting search.")
            return pd.DataFrame()

        final_df, filtered_df, ensp_without_gene_name_df = (
            self._convert_protein_to_gene_names(protein_file_path)
        )
        if final_df.empty:
            logger.error("Failed to convert proteins to gene names.")
            return pd.DataFrame()

        # 将最终结果保存到指定路径
        final_df.to_csv(output_raw_file, index=False)
        filtered_df.to_csv(output_filtered_file, index=False)
        ensp_without_gene_name_df.to_csv(
            output_ensp_without_gene_name_file, index=False
        )
        logger.info(f"Final result saved to {output_raw_file}")
        logger.info(f"Filtered result saved to {output_filtered_file}")
        logger.info(
            f"ENSP without gene name saved to {output_ensp_without_gene_name_file}"
        )

        return final_df


stich_inchikey_target = TCMDataProcessor().search
# 使用示例
if __name__ == "__main__":
    processor = TCMDataProcessor()
    result_df = processor.search(
        inchikey_file="biorange/data/inchikey.csv",
        chemical_file="biorange/data/TCMSP_NGM_STITCH_INCHIKEY.tsv.gz",
        protein_file="biorange/data/9606.protein_chemical.links.transfer.v5.0.tsv.gz",
        output_raw_file="./results/output2/moltarget/stitch_target_raw.csv",
        output_filtered_file="./results/output2/moltarget/stitch_output_target.csv",
        output_ensp_without_gene_name_file="./results/output2/moltarget/ENSP_without_gene_name.csv",
    )
