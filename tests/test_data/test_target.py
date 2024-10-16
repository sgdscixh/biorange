# 测试chembl
# from biorange.target_predict.mol_target import chembl_local
# import pandas as pd

# compound_input = pd.read_csv(
#     "/home/liuyan/projects/package/biorange/biorange/data/inchikey.csv"
# )
# # compound_input = [
# #     {"inchikey": "UCMIRNVEIXFBKS-UHFFFAOYSA-N"},
# #     {"inchikey": "YCIMNLLNPGFGHC-UHFFFAOYSA-N"},
# # ]

# result = chembl_local.chembl_target_prediction(compound_input)
# result.to_csv("chembl_local_result.csv", index=False)
# print(result)


# 测试stitch
from biorange.target_predict.mol_target import stitch_inchikey
import pandas as pd

processor = stitch_inchikey.TCMDataProcessor()
df = pd.read_csv("/home/liuyan/projects/package/biorange/biorange/data/inchikey.csv")
df_list = df["inchikey"].tolist()
result_df = processor.search(inchikeys=df_list)
result_df = result_df  # 去除重复行

result_df.to_csv("./stitch_target.csv", index=False)  # .drop_duplicates()
