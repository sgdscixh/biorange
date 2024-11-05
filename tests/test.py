# from biorange.target_predict.data_processing import admet_filter
# import pandas as pd


# # admet_filter = admet_filter()
# input_file_path = pd.read_csv(
#     "/home/liuyan/projects/package/biorange/biorange/data/first_input_data.csv"
# )
# admet_filtered_ingredients = admet_filter(
#     input_file_path,
#     apply_filter=True,
# )
# print(admet_filtered_ingredients)
# admet_filtered_ingredients.to_csv("./sss.csv")

# # 测试OMIM
# from biorange.target_predict.disease_target.omim import omim_disease_target

# result_df = omim_disease_target("lung cancer")
# print(result_df)
# result_df.to_csv("./omim.csv")

# # 测试ttd
# from biorange.target_predict.disease_target.ttd import ttd_disease_target

# result_df = ttd_disease_target("Osteoarthritis")
# result_df.to_csv("./ttd.csv", index=False)
# print(result_df)


# # 测试stitch
# from biorange.target_predict.mol_target.stitch_inchikey import stich_inchikey_target

# df = pd.read_csv(
#     "/home/liuyan/projects/package/biorange/biorange/data/admet_filtered_ingredients.csv"
# )
# df_list = df["inchikey"].tolist()
# result_df2 = stich_inchikey_target(inchikeys=df_list, combined_score_threshold=200)
# print(result_df2)
# result_df2.to_csv("./stitch.csv", index=False)

# # 测试chembl
# from biorange.target_predict.mol_target.chembl_local import chembl_inchikey_target

# compound_input = ["UCMIRNVEIXFBKS-UHFFFAOYSA-N", "YCIMNLLNPGFGHC-UHFFFAOYSA-N"]
# result = chembl_inchikey_target(compound_input)
# result.to_csv("./chembl.csv")
# print(result)


# #
# from biorange.ppi.generate_type_file import generate_type
# import pandas as pd

# kegg = pd.read_csv("/home/liuyan/projects/package/biorange/biorange/data/kegg_df.csv")
# tar = pd.read_csv("/home/liuyan/projects/package/biorange/biorange/data/target_df.csv")
# node_relationships_df, node_types_df, _ = generate_type(
#     kegg,
#     tar,
# )
# node_relationships_df.to_csv("./node_file.csv", index=False)
# node_types_df.to_csv("./type_file.csv", index=False)

# from biorange.ppi.ppi_final import main
# import os

# gene_names_file = "biorange/data/shared targets of drugs and diseases.csv"
# gene_names = pd.read_csv(gene_names_file)["shared_targets"]

# output_dir = "./results/output2/ppi"
# os.makedirs(output_dir, exist_ok=True)

# # 完整ppi
# results = main(gene_names, output_dir)
# print(results)


# from biorange.ppi.ppi_plot import draw_concentric_layout, draw_custom_layout

# nodes_df = pd.read_csv("/home/liuyan/projects/package/biorange/node_file.csv")
# types_df = pd.read_csv("/home/liuyan/projects/package/biorange/type_file.csv")

# # 使用自定义布局
# draw_custom_layout(
#     nodes_df,
#     types_df,
#     num_rows=7,
#     num_cols=12,
#     output_file="custom_layout",
#     output_dir="./results/output/custom_layout",
#     figsize=(14, 10),
# )

# # 使用同心圆布局
# target_layers = [15, 23, 35, 80]  # 每层的target节点数量
# layer_radii = [
#     0.13,
#     0.35,
#     0.48,
#     0.63,
#     0.78,
#     1.3,
# ]  # 层数为target节点层数+2（compounds节点和pathway节点）

# draw_concentric_layout(
#     nodes_df,
#     types_df,
#     target_layers=target_layers,
#     layer_radii=layer_radii,
#     output_file="concentric_layout",
#     output_dir="./results/output/concentric_layout",
#     figsize=(12, 12),
# )

from biorange.venn.venn_plot import vennplot

TCMSP = ["A", "B", "C", "D", "D", "D"]
Chembl = ["B", "C", "E", "F"]
STITCH = ["C", "D", "E", "G"]

# 使用默认配置
plotter = vennplot(
    [TCMSP, Chembl, STITCH],
    ["TCMSP", "Chembl", "STITCH"],
    "Ingredients_Targets_venn",
    "Ingredients_Targets_venn.png",
)
