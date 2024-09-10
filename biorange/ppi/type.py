import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os


def create_custom_layout(node_types, num_rows, num_cols):
    pos = {}
    y_offsets = {t: -len(node_types[t]) / 2 for t in node_types}

    for i, node in enumerate(node_types["pathway"]):
        pos[node] = (-1.6, i + y_offsets["pathway"])

    for i, node in enumerate(node_types["target"]):
        x = (i % num_cols) * 0.08 - (num_cols * 0.2) / 2
        y = (i // num_cols) * 0.35 + y_offsets["target"]
        pos[node] = (x, y)

    for i, node in enumerate(node_types["counpounds"]):
        pos[node] = (0.3, i + y_offsets["counpounds"])

    centers = {
        t: sum(pos[n][1] for n in node_types[t]) / len(node_types[t])
        for t in node_types
    }
    overall_center = sum(centers.values()) / len(centers)

    for t in node_types:
        offset = overall_center - centers[t]
        for n in node_types[t]:
            pos[n] = (pos[n][0], pos[n][1] + offset)

    return pos


def create_concentric_layout(G, layers, layer_radii):
    pos = {}
    for i, layer in enumerate(layers):
        radius = layer_radii[i]
        angle_step = 2 * np.pi / len(layer)
        for j, node in enumerate(layer):
            angle = j * angle_step
            pos[node] = (radius * np.cos(angle), radius * np.sin(angle))
    return pos


def draw_custom_layout(
    node_file,
    type_file,
    num_rows=7,
    num_cols=12,
    output_file="output",
    figsize=(14, 10),
):
    nodes_df = pd.read_csv(node_file)
    types_df = pd.read_csv(type_file)

    G = nx.Graph()
    for _, row in nodes_df.iterrows():
        G.add_edge(row["node1"], row["node2"])

    type_dict = types_df.set_index("node")["type"].to_dict()
    nx.set_node_attributes(G, type_dict, "type")

    type_color = {"counpounds": "#0D71BF", "target": "#2BB11E", "pathway": "#F56327"}
    type_shape = {"counpounds": "d", "target": "o", "pathway": "*"}

    node_types = {
        t: [n for n, d in G.nodes(data=True) if d["type"] == t] for t in type_color
    }

    pos = create_custom_layout(node_types, num_rows, num_cols)

    plt.figure(figsize=figsize)

    node_sizes = {"counpounds": 200, "target": 300, "pathway": 200}
    font_sizes = {"counpounds": 10, "target": 8, "pathway": 10}

    for node_type, nodes in node_types.items():
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=nodes,
            node_color=type_color[node_type],
            node_shape=type_shape[node_type],
            label=node_type,
            node_size=node_sizes[node_type],
        )
        nx.draw_networkx_labels(
            G, pos, labels={n: n for n in nodes}, font_size=font_sizes[node_type]
        )

    nx.draw_networkx_edges(G, pos, edge_color="grey", width=0.3, alpha=0.45)

    plt.legend(scatterpoints=1, markerscale=0.4, fontsize=8)

    output_dir = "./results/output/type"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.savefig(os.path.join(output_dir, f"{output_file}.pdf"))
    plt.savefig(os.path.join(output_dir, f"{output_file}.png"))
    plt.show()


def draw_concentric_layout(
    node_file,
    type_file,
    target_layers=None,
    layer_radii=None,
    output_file="output",
    figsize=(12, 12),
):
    nodes_df = pd.read_csv(node_file)
    types_df = pd.read_csv(type_file)

    G = nx.Graph()
    for _, row in nodes_df.iterrows():
        G.add_edge(row["node1"], row["node2"])

    type_dict = types_df.set_index("node")["type"].to_dict()
    nx.set_node_attributes(G, type_dict, "type")

    type_color = {"counpounds": "#0D71BF", "target": "#2BB11E", "pathway": "#F56327"}
    type_shape = {"counpounds": "d", "target": "o", "pathway": "*"}

    node_types = {
        t: [n for n, d in G.nodes(data=True) if d["type"] == t] for t in type_color
    }

    if layer_radii is None:
        layer_radii = [0.1, 0.4, 0.55, 0.7, 0.85, 1, 1.6]
    if target_layers is None:
        target_nodes_sorted = sorted(node_types["target"])
        target_layers = [
            target_nodes_sorted[:12],
            target_nodes_sorted[12:30],
            target_nodes_sorted[30:51],
            target_nodes_sorted[51:80],
            target_nodes_sorted[80:130],
        ]
    else:
        target_nodes_sorted = sorted(node_types["target"])
        target_layers = [
            target_nodes_sorted[sum(target_layers[:i]) : sum(target_layers[: i + 1])]
            for i in range(len(target_layers))
        ]

    compounds_nodes = node_types["counpounds"]
    pathway_nodes = node_types["pathway"]
    layers = [compounds_nodes, *target_layers, pathway_nodes]
    pos = create_concentric_layout(G, layers, layer_radii)

    plt.figure(figsize=figsize)

    node_sizes = {"counpounds": 200, "target": 300, "pathway": 200}
    font_sizes = {"counpounds": 10, "target": 8, "pathway": 10}

    for node_type, nodes in node_types.items():
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=nodes,
            node_color=type_color[node_type],
            node_shape=type_shape[node_type],
            label=node_type,
            node_size=node_sizes[node_type],
        )
        nx.draw_networkx_labels(
            G, pos, labels={n: n for n in nodes}, font_size=font_sizes[node_type]
        )

    nx.draw_networkx_edges(G, pos, edge_color="grey", width=0.3, alpha=0.45)

    plt.legend(scatterpoints=1, markerscale=0.4, fontsize=8)

    output_dir = "./results/output/type"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.savefig(os.path.join(output_dir, f"{output_file}.pdf"))
    plt.savefig(os.path.join(output_dir, f"{output_file}.png"))
    plt.show()


if __name__ == "__main__":
    nodes_df_custom = (
        "/home/liuyan/projects/package/biorange/biorange/data/node_file333.csv"
    )
    type_df_custom = (
        "/home/liuyan/projects/package/biorange/biorange/data/type_file.csv"
    )
    nodes_df_concentric = (
        "/home/liuyan/projects/package/biorange/biorange/data/node_file333.csv"
    )
    type_df_concentric = (
        "/home/liuyan/projects/package/biorange/biorange/data/type_file.csv"
    )

    # 使用自定义布局
    draw_custom_layout(
        nodes_df_custom,
        type_df_custom,
        num_rows=7,
        num_cols=12,
        output_file="custom_layout",
        figsize=(14, 10),
    )

    # 使用同心圆布局
    target_layers = [12, 18, 21, 29, 50]  # 每层的节点数量
    layer_radii = [0.1, 0.4, 0.55, 0.7, 0.85, 1, 1.6]
    draw_concentric_layout(
        nodes_df_concentric,
        type_df_concentric,
        target_layers=target_layers,
        layer_radii=layer_radii,
        output_file="concentric_layout",
        figsize=(12, 12),
    )
