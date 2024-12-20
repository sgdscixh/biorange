import pandas as pd
from biorange.utils.package_fileload import get_data_file_path


def chembl_inchikey_target(
    compound_input,
    organism="Homo sapiens",
    confidence_column="confidence80",
    confidence_types=("active", "both"),
    threshold=5,
):

    # 读取大表数据
    chembl_large = get_data_file_path(
        "chembl_25_targets_internal_data_homo_202410.csv.gz"
    )
    chembl_large_table = pd.read_csv(chembl_large)

    # 进行inchikeys筛选
    inchikey_list = compound_input
    chembl_targets = chembl_large_table[
        chembl_large_table["inchikey"].isin(inchikey_list)
    ]

    # 增加筛选条件
    if organism is not None:
        chembl_targets = chembl_targets[chembl_targets["organism"] == organism]

    chembl_targets = chembl_targets[
        (chembl_targets[confidence_column].isin(confidence_types))
        & (chembl_targets["threshold"] >= threshold)
    ]

    return chembl_targets


# 示例调用
if __name__ == "__main__":

    # compound_input = pd.read_csv(
    #     "/home/liuyan/projects/package/biorange/biorange/data/first_input_data.csv"
    # )
    compound_input = ["UCMIRNVEIXFBKS-UHFFFAOYSA-N", "YCIMNLLNPGFGHC-UHFFFAOYSA-N"]

    result = chembl_inchikey_target(compound_input)
    print(result)
