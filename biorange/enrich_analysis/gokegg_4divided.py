import gseapy as gp
import matplotlib.pyplot as plt
from gseapy import barplot, dotplot
import os


def kegg_enrichment_analysis(gene_list, output_dir):
    """
    进行KEGG富集分析并保存结果
    """
    enr_KEGG = gp.enrichr(
        gene_list=gene_list,
        gene_sets=["KEGG_2021_Human"],
        organism="human",
        cutoff=1,
    )
    os.makedirs(output_dir, exist_ok=True)
    enr_KEGG.results.to_csv(os.path.join(output_dir, "enr_KEGG.csv"), index=False)
    return enr_KEGG


def plot_kegg(enr_KEGG, output_dir):
    """
    绘制KEGG富集分析结果的图并保存
    """
    os.makedirs(output_dir, exist_ok=True)

    # 绘制柱状图
    barplot_fig = barplot(
        enr_KEGG.res2d, title="KEGG_2021_Human", figsize=(6, 7), color="darkred"
    )
    barplot_fig.figure.savefig(
        os.path.join(output_dir, "KEGG_barplot.png"), bbox_inches="tight"
    )

    # 绘制点图
    dotplot_fig = dotplot(
        enr_KEGG.res2d,
        title="KEGG_2021_Human",
        cmap="viridis_r",
        size=10,
        figsize=(6, 7),
    )
    dotplot_fig.figure.savefig(
        os.path.join(output_dir, "KEGG_dotplot.png"), bbox_inches="tight"
    )


def go_enrichment_analysis(gene_list, output_dir):
    """
    进行GO富集分析并保存结果
    """
    enr_GO = gp.enrichr(
        gene_list=gene_list,
        gene_sets=[
            "GO_Biological_Process_2021",
            "GO_Molecular_Function_2021",
            "GO_Cellular_Component_2021",
        ],
        organism="Human",
        cutoff=1,
    )
    os.makedirs(output_dir, exist_ok=True)
    enr_GO.results.to_csv(os.path.join(output_dir, "enr_GO.csv"), index=False)
    return enr_GO


def plot_go(enr_GO, output_dir):
    """
    绘制GO富集分析结果的图并保存
    """
    os.makedirs(output_dir, exist_ok=True)

    # 绘制柱状图，交换横纵坐标
    barplot_fig = barplot(
        enr_GO.results,
        column="Adjusted P-value",
        cutoff=1,
        group="Gene_set",
        size=10,
        top_term=5,
        figsize=(7, 5),
        color=["darkred", "darkblue", "green"],
        orientation="horizontal",  # 设置为水平柱状图
    )
    barplot_fig.figure.savefig(
        os.path.join(output_dir, "GO_barplot.png"), bbox_inches="tight"
    )

    # 绘制点图
    dotplot_fig = dotplot(
        enr_GO.results,
        column="Adjusted P-value",
        x="Gene_set",
        cutoff=1,
        size=10,
        top_term=5,
        figsize=(7, 5),
        title="GO",
        xticklabels_rot=45,
        show_ring=False,
        marker="o",
    )
    dotplot_fig.figure.savefig(
        os.path.join(output_dir, "GO_dotplot.png"), bbox_inches="tight"
    )


# 示例数据
gene_list = [
    "MLH1",
    "ECM10",
    "RLI1",
    "SSB1",
    "SSB2",
    "MSH2",
    "STAT3",
    "IL6",
    "S100A9",
    "S100A8",
    "ARG1",
    "NOS2",
]

# 进行KEGG富集分析并绘图
kegg_output_dir = "./KEGG"
enr_KEGG = kegg_enrichment_analysis(gene_list, kegg_output_dir)
plot_kegg(enr_KEGG, kegg_output_dir)

# 进行GO富集分析并绘图
go_output_dir = "./GO"
enr_GO = go_enrichment_analysis(gene_list, go_output_dir)
plot_go(enr_GO, go_output_dir)
