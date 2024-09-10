import pandas as pd
import numpy as np
from plotnine import *


def go_bar_plot(
    enrich_GO, shownNumber=10, color="P-value", bp="red", cc="green", mf="blue"
):
    """
    Create a bar plot for GO enrichment analysis results.
    Args:
        enrich_GO (pd.DataFrame): DataFrame containing GO enrichment results.
        shownNumber (int): Number of top terms to show for each GO category. Default is 10.
        color (str): Column name to use for color intensity. Default is "P-value".
        bp (str): Color for Biological Process. Default is "red".
        cc (str): Color for Cellular Component. Default is "green".
        mf (str): Color for Molecular Function. Default is "blue".
    Returns:
        plotnine.ggplot: A ggplot object representing the GO enrichment bar plot.
    """
    # Select top terms for each GO category
    top_enrich_GO = (
        enrich_GO.groupby("Gene_set")
        .apply(lambda x: x.nsmallest(shownNumber, "P-value"))
        .reset_index(drop=True)
    )

    # Create the plot
    p = (
        ggplot(
            top_enrich_GO,
            aes(
                x="Description",
                y="Count",
                fill="Gene_set",
                color="Gene_set",
                alpha=f'-np.log10(top_enrich_GO["{color}"])',
            ),
        )
        + geom_bar(stat="identity", width=0.5)
        + scale_color_manual(values=[bp, cc, mf])
        + scale_fill_manual(values=[bp, cc, mf])
        + scale_alpha_continuous(name="-log10(P-value)")
        + theme_bw()
        + expand_limits(x=0, y=0)
        + theme(
            panel_border=element_rect(color="grey"),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            legend_box="horizontal",
            legend_position="top",
            axis_text_x=element_text(
                color="black", size=16, angle=80, hjust=1, family="serif"
            ),
            axis_text_y=element_text(color="black", size=16, family="serif"),
            legend_title=element_text(color="black", size=18, family="serif"),
            legend_text=element_text(color="black", size=16, family="serif"),
            figure_size=(20, 10),
        )
        + guides(alpha=guide_legend(order=2))
        + labs(x="GO Term", y="Gene Count")
    )
    return p


def kegg_bar_plot(
    dat, shownNumber=10, color="P-value", kegg_plot_path=".", fc="#9966CC"
):
    # 选择展示的行数并直接从 'Overlap' 列提取 Count
    enrich_kegg = dat.head(min(shownNumber, len(dat))).copy()
    enrich_kegg["Count"] = enrich_kegg["Overlap"].str.split("/").str[0].astype(int)

    # 使用 'Term' 列作为 'Description'
    enrich_kegg["Description"] = pd.Categorical(
        enrich_kegg["Term"], categories=enrich_kegg["Term"], ordered=True
    )

    # 创建绘图
    p = (
        ggplot(
            enrich_kegg,
            aes(x="Description", y=f'-np.log10(enrich_kegg["{color}"])', alpha="Count"),
        )
        + geom_bar(stat="identity", color=fc, fill=fc, width=0.5)
        + coord_flip()
        + expand_limits(x=0, y=0)
        + theme_classic(base_size=12)
        + xlab(" ")
        + ylab(f"-log10({color})")
        + theme(
            panel_background=element_blank(),
            panel_grid_minor=element_blank(),
            panel_grid_major=element_blank(),
            axis_text_x=element_text(
                face="bold", color="black", size=26, family="serif"
            ),
            axis_text_y=element_text(
                face="bold", color="black", size=26, family="serif"
            ),
            axis_title_x=element_text(
                face="bold", color="black", size=28, family="serif"
            ),
            axis_title_y=element_text(
                face="bold", color="black", size=28, family="serif"
            ),
            legend_title=element_text(color="black", size=26, family="serif"),
            legend_text=element_text(color="black", size=22, family="serif"),
            figure_size=(10, 10),
        )
    )

    return p
