import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.cm as cm
from io import StringIO


def fetch_ppi_data(gene_names, species_id=9606):
    url = "https://string-db.org/api/tsv/network"
    params = {
        "identifiers": "\r".join(gene_names),
        "species": species_id,
        "caller_identity": "my_script",
    }

    response = requests.post(url, data=params)

    if response.status_code != 200:
        print(f"Error: {response.status_code}")
        print(response.text)
        return None

    data = response.text
    return data


def parse_interaction_data(data):
    try:
        interaction_data = pd.read_csv(StringIO(data), sep="\t")
        if (
            "preferredName_A" in interaction_data.columns
            and "preferredName_B" in interaction_data.columns
        ):
            nodes = interaction_data[["preferredName_A", "preferredName_B"]]
            return nodes
        else:
            print("Expected columns not found in the response.")
            return None
    except pd.errors.EmptyDataError:
        print("No data received from STRING API.")
        return None


def save_interaction_data(nodes, output_dir):
    output_csv_path = os.path.join(output_dir, "protein_interactions.csv")
    nodes.to_csv(output_csv_path, index=False, header=["node1", "node2"])
    return output_csv_path


def plot_ppi_network(nodes, output_dir, edge_width=0.1, dpi=900):
    G = nx.from_pandas_edgelist(nodes, "preferredName_A", "preferredName_B")
    plt.figure(figsize=(10, 10))
    pos = nx.kamada_kawai_layout(G)
    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=30,
        node_color="#4562E0",
        font_size=9,
        edge_color="gray",
        width=edge_width,
    )
    plt.title("Protein-Protein Interaction Network")
    png_path = os.path.join(output_dir, "ppi_network.png")
    pdf_path = os.path.join(output_dir, "ppi_network.pdf")
    plt.savefig(png_path, dpi=dpi)
    plt.savefig(pdf_path)
    plt.show()
    return png_path, pdf_path


def calculate_node_degrees(nodes, output_dir):
    G = nx.Graph()
    G.add_edges_from(zip(nodes["preferredName_A"], nodes["preferredName_B"]))
    degree_dict = dict(G.degree(G.nodes()))
    df = pd.DataFrame(list(degree_dict.items()), columns=["node", "degree"])
    df_sorted = df.sort_values(by="degree", ascending=False)
    degree_csv_path = os.path.join(output_dir, "string_node_degree.csv")
    df_sorted.to_csv(degree_csv_path, index=False)
    return df_sorted, degree_csv_path


def plot_core_targets(degree_df, nodes, output_dir):
    degree_dict = degree_df.set_index("node")["degree"].to_dict()
    edges = [
        (row["preferredName_A"], row["preferredName_B"]) for _, row in nodes.iterrows()
    ]

    G = nx.Graph()
    G.add_nodes_from(degree_dict.keys())
    G.add_edges_from(edges)

    top_nodes = [
        node
        for node, degree in sorted(
            degree_dict.items(), key=lambda x: x[1], reverse=True
        )[:274]
    ]
    subgraph = G.subgraph(top_nodes)

    sorted_nodes = sorted(
        {node: degree_dict[node] for node in top_nodes}.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    layers = [1, 6, 13, 20, 30, 40, 44, 55, 65]
    radii = [0, 0.8, 1.5, 2.15, 2.8, 3.4, 3.9, 4.4, 4.8]

    pos = {}
    layer_start = 0
    for layer, radius in zip(layers, radii):
        layer_nodes = sorted_nodes[layer_start : layer_start + layer]
        theta = np.linspace(0, 2 * np.pi, len(layer_nodes), endpoint=False)
        for (node, degree), angle in zip(layer_nodes, theta):
            pos[node] = (radius * np.cos(angle), radius * np.sin(angle))
        layer_start += layer

    degrees = [degree_dict[node] for node in subgraph.nodes()]
    vmin = np.percentile(degrees, 0.005)
    vmax = np.percentile(degrees, 100)
    norm = plt.Normalize(vmin, vmax)
    cmap = cm.plasma
    node_colors = [cmap(norm(degree_dict[node])) for node in subgraph.nodes()]

    min_size = 100
    max_size = 1000
    node_sizes = [
        min_size
        + (max_size - min_size)
        * (degree_dict[node] - min(degrees))
        / (max(degrees) - min(degrees))
        for node in subgraph.nodes()
    ]

    min_font_size = 4
    max_font_size = 10
    font_sizes = {
        node: min_font_size
        + (max_font_size - min_font_size)
        * (degree_dict[node] - min(degrees))
        / (max(degrees) - min(degrees))
        for node in subgraph.nodes()
    }

    plt.figure(figsize=(6, 6))
    nx.draw(
        subgraph,
        pos,
        with_labels=False,
        node_size=node_sizes,
        node_color=node_colors,
        edge_color="gray",
    )
    nx.draw_networkx_edges(subgraph, pos, edge_color="gray")

    for node, (x, y) in pos.items():
        plt.text(x, y, s=node, fontsize=font_sizes[node], ha="center", va="center")

    plt.title("PPI network top 100")
    pdf_path = os.path.join(output_dir, "PPI_network_degree.pdf")
    png_path = os.path.join(output_dir, "PPI_network_degree.png")
    plt.savefig(pdf_path)
    plt.savefig(png_path)
    plt.show()
    return png_path, pdf_path


def main(gene_names, output_dir):
    # 获取PPI数据
    data = fetch_ppi_data(gene_names)
    if data is None:
        return

    # 解析PPI数据
    nodes = parse_interaction_data(data)
    if nodes is None:
        return

    # 保存节点关系数据
    interaction_csv_path = save_interaction_data(nodes, output_dir)
    if interaction_csv_path is None:
        return

    # 绘制PPI网络图
    ppi_png_path, ppi_pdf_path = plot_ppi_network(nodes, output_dir)

    # 计算节点度数并保存
    degree_df, degree_csv_path = calculate_node_degrees(nodes, output_dir)

    # 绘制核心靶点图
    core_png_path, core_pdf_path = plot_core_targets(degree_df, nodes, output_dir)

    return {
        "interaction_csv": interaction_csv_path,
        "ppi_network_png": ppi_png_path,
        "ppi_network_pdf": ppi_pdf_path,
        "degree_csv": degree_csv_path,
        "core_targets_png": core_png_path,
        "core_targets_pdf": core_pdf_path,
    }


if __name__ == "__main__":
    # 输入示例
    gene_names_file = "biorange/data/shared targets of drugs and diseases.csv"
    gene_names = pd.read_csv(gene_names_file)["shared_targets"]

    output_dir = "./results/output2/ppi"
    os.makedirs(output_dir, exist_ok=True)

    # 完整ppi
    results = main(gene_names, output_dir)
    print(results)
