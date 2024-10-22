from biorange.ppi import ppi_final
import pandas as pd


class NetworkTypeProcessor:
    def __init__(self, species_id=9606):
        self.species_id = species_id

    def _extract_kegg_attributes(self, kegg_df):
        kegg_df["gene_name"] = kegg_df["gene_name"].str.split(";")
        expanded_kegg_df = kegg_df.explode("gene_name")

        gene_description_df = expanded_kegg_df[["gene_name", "Description"]]
        return gene_description_df

    def _fetch_ppi_data(self, gene_names):
        ppi_data = ppi_final.fetch_ppi_data(gene_names)
        if ppi_data is None:
            return pd.DataFrame()

        interaction_nodes = ppi_final.parse_interaction_data(ppi_data)
        if interaction_nodes is None:
            return pd.DataFrame()

        # Filter interaction nodes to only include those present in gene_names
        valid_gene_names = set(gene_names)
        filtered_nodes = interaction_nodes[
            interaction_nodes["preferredName_A"].isin(valid_gene_names)
            & interaction_nodes["preferredName_B"].isin(valid_gene_names)
        ]

        filtered_nodes.to_csv("./inter.csv", index=False)
        return filtered_nodes

    def _associate_ingredients(self, gene_description_df, targets_df):
        merged_df = gene_description_df.merge(
            targets_df[["gene_name", "compound_name"]],
            on="gene_name",
            how="left",
        )
        gene_ingredient_df = merged_df[["gene_name", "compound_name"]].dropna()

        return gene_ingredient_df

    def _create_output_files(
        self, gene_description_df, gene_ingredient_df, interaction_nodes
    ):
        node_relationships_df = pd.concat(
            [
                gene_description_df[["gene_name", "Description"]].rename(
                    columns={"gene_name": "node1", "Description": "node2"}
                ),
                gene_ingredient_df[["gene_name", "compound_name"]].rename(
                    columns={"gene_name": "node1", "compound_name": "node2"}
                ),
                interaction_nodes[["preferredName_A", "preferredName_B"]].rename(
                    columns={"preferredName_A": "node1", "preferredName_B": "node2"}
                ),
            ],
            ignore_index=True,
        )

        node_types_df = pd.concat(
            [
                gene_description_df[["Description"]]
                .rename(columns={"Description": "node"})
                .assign(type="pathway"),
                gene_description_df[["gene_name"]]
                .rename(columns={"gene_name": "node"})
                .assign(type="target"),
                gene_ingredient_df[["compound_name"]]
                .rename(columns={"compound_name": "node"})
                .assign(type="compound"),
            ],
            ignore_index=True,
        ).drop_duplicates()

        return node_relationships_df, node_types_df

    def process_from_file(self, kegg_file_path, targets_total_file):
        kegg_df = pd.read_csv(kegg_file_path)
        targets_total_df = pd.read_csv(targets_total_file)
        return self._process(kegg_df, targets_total_df)

    def process_from_dataframe(self, kegg_df, targets_total_df):
        return self._process(kegg_df, targets_total_df)

    def _process(self, kegg_df, targets_total_df):
        gene_description_df = self._extract_kegg_attributes(kegg_df)
        interaction_nodes = self._fetch_ppi_data(
            gene_description_df["gene_name"].unique()
        )
        gene_ingredient_df = self._associate_ingredients(
            gene_description_df, targets_total_df
        )
        node_relationships_df, node_types_df = self._create_output_files(
            gene_description_df, gene_ingredient_df, interaction_nodes
        )
        return node_relationships_df, node_types_df, interaction_nodes


generate_type = NetworkTypeProcessor().process_from_dataframe

if __name__ == "__main__":
    processor = NetworkTypeProcessor()
    node_relationships_df, node_types_df, _ = processor.process_from_file(
        "/home/liuyan/projects/package/biorange/biorange/data/kegg_df.csv",
        "/home/liuyan/projects/package/biorange/biorange/data/target_df.csv",
    )
    node_relationships_df.to_csv("node_file88.csv", index=False)
    node_types_df.to_csv("type_file88.csv", index=False)
