import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from biorange.venn.venn_config import VennPlotConfig


class VennPlotter:
    def __init__(self, config: VennPlotConfig = None):
        self.config = config or VennPlotConfig()
        plt.rcParams["font.family"] = [self.config.font_family]

    def plot_venn(self, groups, labels, title, filename=None, **kwargs):
        self.config.update(**kwargs)  # Update config with any additional kwargs
        fig, ax = plt.subplots(figsize=self.config.figsize, dpi=self.config.dpi)

        if len(groups) == 2:
            vee = venn2(
                [set(groups[0]), set(groups[1])],
                set_labels=labels,
                set_colors=self.config.set_colors[:2],
                alpha=self.config.alpha,
                ax=ax,
            )
            venn2_circles(
                [set(groups[0]), set(groups[1])],
                linestyle=self.config.linestyle,
                linewidth=self.config.linewidth,
                color=self.config.circle_color,
                ax=ax,
            )
        elif len(groups) == 3:
            vee = venn3(
                [set(groups[0]), set(groups[1]), set(groups[2])],
                set_labels=labels,
                set_colors=self.config.set_colors,
                alpha=self.config.alpha,
                ax=ax,
            )
            venn3_circles(
                [set(groups[0]), set(groups[1]), set(groups[2])],
                linestyle=self.config.linestyle,
                linewidth=self.config.linewidth,
                color=self.config.circle_color,
                ax=ax,
            )
        else:
            raise ValueError("Only 2 or 3 groups are supported.")

        # Set labels' font properties
        for text in vee.set_labels or []:
            if text:
                text.set_fontsize(self.config.label_fontsize)
                text.set_fontweight(self.config.label_fontweight)
        for text in vee.subset_labels or []:
            if text:
                text.set_fontsize(self.config.label_fontsize)
                text.set_fontweight(self.config.label_fontweight)

        # Set plot title with custom font properties
        plt.title(
            title,
            fontsize=self.config.title_fontsize,
            fontweight=self.config.title_fontweight,
            color=self.config.title_color,
            style=self.config.title_style,
        )

        if filename:
            output_dir = os.path.join("results", "venn")
            os.makedirs(output_dir, exist_ok=True)
            filepath = os.path.join(output_dir, f"venn-{filename}")
            print(f"Saving plot to {filepath}")
            plt.savefig(filepath)
            self.intersection(groups, labels, filename)
        self.intersection(groups, labels)
        plt.show()

    def intersection(self, groups, labels, filename=None):
        sets = [set(group) for group in groups]
        intersections = {}

        if len(sets) == 2:
            intersections[f"{labels[0]}&{labels[1]}"] = sets[0] & sets[1]
        elif len(sets) == 3:
            intersections[f"{labels[0]}&{labels[1]}"] = sets[0] & sets[1]
            intersections[f"{labels[0]}&{labels[2]}"] = sets[0] & sets[2]
            intersections[f"{labels[1]}&{labels[2]}"] = sets[1] & sets[2]
            intersections[f"{labels[0]}&{labels[1]}&{labels[2]}"] = (
                sets[0] & sets[1] & sets[2]
            )

        # 将交集转换为 pandas DataFrame，交集名称作为列名
        df = pd.DataFrame(
            dict(
                [(key, pd.Series(list(value))) for key, value in intersections.items()]
            )
        )

        if filename:
            output_dir = os.path.join("results", "venn")
            filepath = os.path.join(
                output_dir, f"venn-data-{os.path.splitext(filename)[0]}.csv"
            )
            df.to_csv(filepath, index=False)
            print(f"Intersection data saved to {filepath}")
        else:
            print(df)


vennplot = VennPlotter().plot_venn

# Usage Example
if __name__ == "__main__":
    # 示例数据
    TCMSP = ["A", "B", "C", "D", "D", "D"]
    Chembl = ["B", "C", "E", "F"]
    STITCH = ["C", "D", "E", "G"]

    # 使用默认配置
    plotter = VennPlotter()
    vennplot(
        [TCMSP, Chembl, STITCH],
        ["TCMSP", "Chembl", "STITCH"],
        "Ingredients_Targets_venn",
        "Ingredients_Targets_venn.png",
    )

    Genecards = ["H", "I", "J", "K"]
    OMIM = ["I", "J", "L", "M"]
    TTD = ["J", "K", "L", "N"]
    plotter.plot_venn(
        [Genecards, OMIM, TTD],
        ["Genecards", "OMIM", "TTD"],
        "Disease_Targets_venn",
        "Disease_Targets_venn.png",
    )

    Ingredient_Targets = ["O", "P", "Q", "R"]
    Disease_Targets = ["P", "Q", "S", "T"]
    plotter.plot_venn(
        [Ingredient_Targets, Disease_Targets],
        ["Ingredient_Targets", "Disease_Targets"],
        "Ingredient_Disease_Targets_veen",
        "Ingredient_Disease_Targets_veen.png",
    )

    # 使用自定义配置
    custom_config = VennPlotConfig(
        figsize=(8, 6), dpi=150, title_fontsize=16, label_fontsize=12
    )
    custom_plotter = VennPlotter(config=custom_config)
    custom_plotter.plot_venn(
        [TCMSP, Chembl],
        ["TCMSP", "Chembl"],
        "Custom Venn Diagram",
        "test_custom_venn.pdf",
    )

    # 动态更新配置
    plotter.plot_venn(
        [TCMSP, Chembl, STITCH],
        ["TCMSP", "Chembl", "STITCH"],
        "Updated Venn Diagram",
        "test_updated_venn.pdf",
        title_fontsize=18,
        label_fontsize=14,
    )
