from biorange.target_predict import ChEMBLTargetScraper


client = ChEMBLTargetScraper()
smiles = "CCO"
df = client.search_smiles(smiles)
df.to_csv("results/output.csv", index=False)

from biorange.venn.venn_plot import VennPlotter

# 示例数据
TCMSP = ["A", "B", "C", "D"]
Chembl = ["B", "C", "E", "F"]
STITCH = ["C", "D", "E", "G"]
# 使用默认配置
plotter = VennPlotter()
plotter.plot_venn(
    [TCMSP, Chembl, STITCH],
    ["TCMSP", "Chembl", "STITCH"],
    "Ingredients_Targets_venn",
    "Ingredients_Targets_venn.png",
)
