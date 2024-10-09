from biorange.ppi import ppi_final
import pandas as pd


class NetworkTypeProcessor:
    def __init__(self, species_id=9606):
        self.species_id = species_id

    def _extract_kegg_attributes(self, kegg_df):
        kegg_df["gene_name"] = kegg_df["gene_name"].str.split(";")
        expanded_kegg_df = kegg_df.explode("gene_name")

        gene_term_df = expanded_kegg_df[["gene_name", "term"]]

        return gene_term_df

    def _fetch_ppi_data(self, gene_names):
        ppi_data = ppi_final.fetch_ppi_data(gene_names)
        if ppi_data is None:
            return
        interaction_nodes = ppi_final.parse_interaction_data(ppi_data)
        if interaction_nodes is None:
            return
        return interaction_nodes

    def _associate_ingredients(self, gene_term_df, targets_df):
        merged_df = gene_term_df.merge(
            targets_df[["gene_name", "ingredient_name"]],
            on="gene_name",
            how="left",
        )
        gene_ingredient_df = merged_df[["gene_name", "ingredient_name"]].dropna()

        return gene_ingredient_df

    def _create_output_files(self, gene_term_df, gene_ingredient_df, interaction_nodes):

        node_relationships_df = pd.concat(
            [
                gene_term_df[["gene_name", "term"]].rename(
                    columns={"gene_name": "node1", "term": "node2"}
                ),
                gene_ingredient_df[["gene_name", "ingredient_name"]].rename(
                    columns={"gene_name": "node1", "ingredient_name": "node2"}
                ),
                interaction_nodes[["preferredName_A", "preferredName_B"]].rename(
                    columns={"preferredName_A": "node1", "preferredName_B": "node2"}
                ),
            ],
            ignore_index=True,
        )

        node_types_df = pd.concat(
            [
                gene_term_df[["term"]]
                .rename(columns={"term": "node"})
                .assign(type="pathway"),
                gene_term_df[["gene_name"]]
                .rename(columns={"gene_name": "node"})
                .assign(type="target"),
                gene_ingredient_df[["ingredient_name"]]
                .rename(columns={"ingredient_name": "node"})
                .assign(type="compound"),
            ],
            ignore_index=True,
        ).drop_duplicates()

        return node_relationships_df, node_types_df

    # 从文件中处理
    def process_from_file(self, kegg_file_path, targets_total_file):
        kegg_df = pd.read_csv(kegg_file_path)
        targets_total_df = pd.read_csv(targets_total_file)
        return self._process(kegg_df, targets_total_df)

    # 从DataFrame中处理，你喜欢从文件处理的话，以后让AI多生成一个纯df的方法，真的太影响了，后面是需要拓展的 从文件处理拓展不了 写死了
    def process_from_dataframe(self, kegg_df, targets_total_df):
        return self._process(kegg_df, targets_total_df)

    def _process(self, kegg_df, targets_total_df):
        gene_term_df = self._extract_kegg_attributes(kegg_df)
        interaction_nodes = self._fetch_ppi_data(gene_term_df["gene_name"].unique())
        gene_ingredient_df = self._associate_ingredients(gene_term_df, targets_total_df)
        node_relationships_df, node_types_df = self._create_output_files(
            gene_term_df, gene_ingredient_df, interaction_nodes
        )
        return node_relationships_df, node_types_df


generate_type = NetworkTypeProcessor().process_from_dataframe


if __name__ == "__main__":
    # 示例调用
    processor = NetworkTypeProcessor()
    node_relationships_df, node_types_df = processor.process_from_file(
        "/home/liuyan/projects/package/biorange/notebooks/kegg_test.csv",
        "/home/liuyan/projects/package/biorange/notebooks/targets_total_file.csv",
    )
    # node_relationships_df.to_csv("node_file33.csv", index=False)
    # node_types_df.to_csv("type_file33.csv", index=False)
