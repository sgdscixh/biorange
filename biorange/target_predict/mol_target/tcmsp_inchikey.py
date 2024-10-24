import pandas as pd
from biorange.logger import get_logger
from biorange.utils.package_fileload import get_data_file_path

logger = get_logger(__name__)


class TCMSPTargetScraper:
    def __init__(self, molecules_csv="TCMSP_mol.csv", targets_csv="TCMSP_tar.csv"):
        self.molecules_df = pd.read_csv(get_data_file_path(molecules_csv))
        self.targets_df = pd.read_csv(get_data_file_path(targets_csv))
        self.merged_df = self._merge_dataframes()

    def _merge_dataframes(self):
        merged_df = self.molecules_df.merge(
            self.targets_df, left_on="molecule_ID", right_on="molecule_ID", how="left"
        )
        logger.debug(f"Merged dataframe shape: {merged_df.shape}")
        return merged_df

    def search_inchikeys(
        self,
        inchikeys,
        internal_data_file: str = get_data_file_path(
            "TCM_NGM_inchike_isosmile_11294.csv"
        ),
    ):
        logger.info(f"Searching for InChIKeys: {inchikeys}")

        # 筛选输入的inchikeys
        filtered_df = self.merged_df[self.merged_df["inchikey"].isin(inchikeys)].copy()
        filtered_df["source"] = "TCMSP"

        if filtered_df.empty:
            logger.warning(f"No matches found for InChIKeys: {inchikeys}")
            # 如果没有匹配的inchikeys
            results_df = pd.DataFrame(
                {
                    "inchikey": inchikeys,
                    "gene_name": [None] * len(inchikeys),
                    "source": ["TCMSP"] * len(inchikeys),
                }
            )
        else:
            logger.info(f"Matches found for InChIKeys: {inchikeys}")
            # 选择需要的列并重命名
            results_df = filtered_df[["inchikey", "Gene Names", "source"]].rename(
                columns={"Gene Names": "gene_name"}
            )
            results_df = results_df.dropna(subset=["gene_name"])

        # 加载内置数据
        internal_data_df = pd.read_csv(
            internal_data_file  # , sep="\t"
        )  # 根据实际文件格式调整
        # 合并内置数据
        merged_df = pd.merge(results_df, internal_data_df, on="inchikey", how="left")
        return merged_df


tcmsp_inchikey_target = TCMSPTargetScraper().search_inchikeys

# 示例使用
if __name__ == "__main__":
    searcher = TCMSPTargetScraper()

    # 从CSV文件读取inchikeys
    inchikeys_df = pd.read_csv(
        "/home/liuyan/projects/package/biorange/biorange/data/inchikey.csv"
    )
    inchikeys_list = inchikeys_df["inchikey"].tolist()

    result = searcher.search_inchikeys(inchikeys_list)
    result.to_csv("./tcmsp_output_target.csv", index=False)
