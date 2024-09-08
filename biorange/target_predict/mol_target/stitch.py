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

    def _merge_smiles_to_inchikey(
        self, smiles_file: str, ngm_file: str
    ) -> pd.DataFrame:
        """将 SMILES 转换为 InChIKey 并输出中间结果。"""
        df_smiles = self._read_csv(smiles_file)
        df_ngm = self._read_csv(ngm_file)
        merged_df = self._merge_files(
            df_smiles, df_ngm[["smiles", "inchikey"]], on="smiles"
        )
        inchikey_file_path = self._save_dataframe(
            merged_df, "TCM_ingredients_inchikey.csv"
        )
        return merged_df, inchikey_file_path

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
        protein_file_path = self._save_dataframe(result_df, "TCM_stitch_targets.csv")
        return result_df, protein_file_path

    def _convert_protein_to_gene_names(self, protein_file: str) -> pd.DataFrame:
        """将 protein ID 转换为基因名。"""
        df = self._read_csv(protein_file)
        if df.empty:
            logger.error(f"{protein_file} is empty. Skipping gene name conversion.")
            return pd.DataFrame()

        df["ENSP"] = df["protein"].str.split(".").str[1]
        df = df.drop_duplicates(subset=["ENSP"])
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
        final_df = final_df[["smiles", "gene_name", "source"]]
        final_file_path = self._save_dataframe(final_df, "STITCH_output_targets.csv")
        return final_df

    def search(
        self, smiles_file: str, ngm_file: str, chemical_file: str, protein_file: str
    ) -> pd.DataFrame:
        """主接口：根据 SMILES 查找对应的基因名。"""
        inchikey_df, inchikey_file_path = self._merge_smiles_to_inchikey(
            smiles_file, ngm_file
        )
        if inchikey_df.empty:
            logger.error("Failed to generate InChIKey data. Aborting search.")
            return pd.DataFrame()

        chemical_df, chemical_file_path = self._merge_with_large_chemical_data(
            inchikey_file_path, chemical_file
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

        final_df = self._convert_protein_to_gene_names(protein_file_path)
        if final_df.empty:
            logger.error("Failed to convert proteins to gene names.")
            return pd.DataFrame()

        return final_df


# 使用示例
if __name__ == "__main__":
    processor = TCMDataProcessor()
    result_df = processor.search(
        "/home/liuyan/projects/netparam/biorange/data/smiles.csv",
        "/home/liuyan/projects/netparam/biorange/data/中药NGM_combined_ingredients_inchikey_smiles.csv",
        "biorange/data/TCMSP_NGM_STITCH_INCHIKEY.tsv.gz",
        "/home/liuyan/projects/netparam/biorange/data/9606.protein_chemical.links.transfer.v5.0.tsv.gz",
    )
    result_df.to_csv("results/stich_output.csv", index=False)
