import pandas as pd
from pathlib import Path

# 内置数据在python中主要是相对位置问题  之前写了一个获取内置数据函数，读取这个包data下指定名字的数据
from biorange.utils.package_fileload import get_data_file_path


class ADMETFilter:
    def __init__(
        self,
        lipinski_threshold=4,
        qed_threshold=0.5,
        bioavailability_threshold=0.3,
    ):
        self.lipinski_threshold = lipinski_threshold
        self.qed_threshold = qed_threshold
        self.bioavailability_threshold = bioavailability_threshold
        self.admet_file_path = get_data_file_path("TCM_NGM__ADMET.csv")
        self.admet_data = self._load_admet_data()

    def _load_admet_data(self):
        # 检查 ADMET 数据文件是否存在
        if not self.admet_file_path.is_file():
            raise FileNotFoundError(f"未找到 ADMET 数据文件：{self.admet_file_path}")

        # 读取 ADMET 数据
        admet_data = pd.read_csv(self.admet_file_path)
        required_columns = ["inchikey", "Lipinski", "QED", "Bioavailability_Ma"]
        missing_columns = [
            col for col in required_columns if col not in admet_data.columns
        ]
        if missing_columns:
            raise ValueError(f"ADMET 数据缺少必要的列：{', '.join(missing_columns)}")

        return admet_data

    def process_file(self, input_file_path, apply_filter=True):
        # 从文件读取输入数据
        input_data = pd.read_csv(input_file_path)
        return self._process_data(input_data, apply_filter)

    def process_dataframe(self, input_dataframe, apply_filter=True):
        # 处理 DataFrame 格式的输入数据
        return self._process_data(input_dataframe, apply_filter)

    def _process_data(self, data, apply_filter):
        # 检查输入数据是否为 DataFrame
        if not isinstance(data, pd.DataFrame):
            raise TypeError("输入数据必须是一个 pandas DataFrame。")

        # 检查输入数据中是否包含 'inchikey' 列
        if "inchikey" not in data.columns:
            raise ValueError("输入数据必须包含 'inchikey' 列。")

        # 合并输入数据与 ADMET 数据
        merged_data = pd.merge(data, self.admet_data, on="inchikey", how="left")

        # 检查合并结果是否为空
        if merged_data.empty:
            raise ValueError(
                "合并后的数据为空，请检查输入数据中的 'inchikey' 列是否有效。"
            )

        if apply_filter:
            # 确保用于过滤的列存在
            filter_columns = ["Lipinski", "QED", "Bioavailability_Ma"]
            missing_filter_columns = [
                col for col in filter_columns if col not in merged_data.columns
            ]
            if missing_filter_columns:
                raise ValueError(
                    f"缺少用于过滤的列：{', '.join(missing_filter_columns)}"
                )

            # 应用过滤条件
            filtered_data = merged_data[
                (merged_data["Lipinski"] == self.lipinski_threshold)
                & (merged_data["QED"] > self.qed_threshold)
                & (merged_data["Bioavailability_Ma"] > self.bioavailability_threshold)
            ]
        else:
            filtered_data = merged_data

        # 选择需要保留的列
        columns_to_keep = [
            "inchikey",
            "smiles",
            "Name",
            "molecular_weight",
            "logP",
            "hydrogen_bond_acceptors",
            "hydrogen_bond_donors",
            "Lipinski",
            "QED",
            "Bioavailability_Ma",
        ]

        # 检查需要保留的列是否存在
        existing_columns = [
            col for col in columns_to_keep if col in filtered_data.columns
        ]
        if not existing_columns:
            raise ValueError("过滤后的数据中没有需要保留的列。")

        return filtered_data[existing_columns]


admet_filter = ADMETFilter().process_dataframe
if __name__ == "__main__":
    # 示例
    input_file_path = (
        "/home/liuyan/projects/package/biorange/biorange/data/first_input_data.csv"
    )

    admet_filter = ADMETFilter()
    admet_filtered_ingredients = admet_filter.process_file(
        input_file_path, apply_filter=True
    )
    print(admet_filtered_ingredients)
