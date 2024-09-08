# 计算节点的代码待传入
# 成功跑出网络图和节点表
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os

# 输入的基因名列表
# gene_names = [
#     "TP53", "AKT1", "EGFR", "GAPDH", "MYC", "CTNNB1", "ESR1", "PTEN", "STAT3", "JUN"
# ]
gene_names = pd.read_csv(
    "/home/liuyan/projects/netparam/biorange/data/PPI/shared targets of drugs and diseases.csv"
)["shared_targets"]
# API URL
url = "https://string-db.org/api/tsv/network"

# 构建POST请求的数据
params = {
    "identifiers": "\r".join(gene_names),
    "species": 9606,  # 9606 是人类的物种ID
    "caller_identity": "my_script",
}

# 发送POST请求
response = requests.post(url, data=params)

# 检查请求是否成功
if response.status_code != 200:
    print(f"Error: {response.status_code}")
    print(response.text)
    exit()

# 打印响应内容以进行调试
print("Response received from STRING API:")
print(response.text[:500])  # 只打印前500个字符以检查响应内容

# 解析响应数据
data = response.text

# 将响应数据转换为DataFrame
from io import StringIO

try:
    interaction_data = pd.read_csv(StringIO(data), sep="\t")
except pd.errors.EmptyDataError:
    print("No data received from STRING API.")
    exit()

# 检查数据框内容
print("DataFrame head:")
print(interaction_data.head())

# 提取节点关系
if (
    "preferredName_A" in interaction_data.columns
    and "preferredName_B" in interaction_data.columns
):
    nodes = interaction_data[["preferredName_A", "preferredName_B"]]
else:
    print("Expected columns not found in the response.")
    exit()

# 获取当前工作目录
current_dir = os.getcwd()
print(f"Current working directory: {current_dir}")

# 检查并创建输出目录
output_dir = os.path.join(current_dir, "output")
os.makedirs(output_dir, exist_ok=True)

# 保存节点关系到CSV文件
output_csv_path = os.path.join(output_dir, "protein_interactions.csv")
try:
    nodes.to_csv(output_csv_path, index=False, header=["node1", "node2"])
    print(f"Protein-protein interactions have been saved to {output_csv_path}")
except PermissionError:
    print(f"Permission denied: Unable to save file to {output_csv_path}")
    exit()

# 再次读取CSV文件以验证其内容
try:
    read_data = pd.read_csv(output_csv_path)
    print("CSV file content:")
    print(read_data.head())
except Exception as e:
    print(f"Error reading the CSV file: {e}")

# 构建PPI网络图
G = nx.from_pandas_edgelist(nodes, "preferredName_A", "preferredName_B")

# 绘制PPI网络图
plt.figure(figsize=(6, 6))
pos = nx.spring_layout(G)
nx.draw(
    G,
    pos,
    with_labels=True,
    node_size=3000,
    node_color="skyblue",
    font_size=10,
    font_weight="bold",
    edge_color="gray",
)
plt.title("Protein-Protein Interaction Network")
plt.savefig(os.path.join(output_dir, "ppi_network.png"))
plt.show()

print(
    f"PPI network graph has been saved to {os.path.join(output_dir, 'ppi_network.png')}"
)


############################


# 计算degree值
# 输入节点表
gf = pd.read_csv("output/protein_interactions.csv")

# 创建图
G = nx.Graph()
G.add_edges_from(zip(gf.iloc[:, 0], gf.iloc[:, 1]))

# 计算每个节点的度数
degree_dict = dict(G.degree(G.nodes()))

# 将度数添加到 DataFrame
df = pd.DataFrame(list(degree_dict.items()), columns=["node", "degree"])

# 根据度数排序
df_sorted2 = df.sort_values(by="degree", ascending=False)
df_sorted2.to_csv("string_node_degree.csv", index=False)

# 读取保存的 CSV 文件
df_read = pd.read_csv("string_node_degree.csv")

# 打印排序后的节点和度数
print(df_read)

# # 设置绘图布局为圆形
# pos = nx.circular_layout(G, scale=0.3)
# # 图有点丑，可以暂时不展示


########################################
# 根据degree文件和节点文件绘制核心靶点图
import numpy as np
import matplotlib.cm as cm

# 读取节点Degree值表
# degree_df = pd.read_csv('string_node_degree.csv')  # 假设文件名为degree_test.csv

degree_df = df_read  # 假设文件名为degree_test.csv
degree_dict = degree_df.set_index("node")["degree"].to_dict()

# 读取边信息表
edges_df = pd.read_csv("output/protein_interactions.csv")  # 假设文件名为edges_test.csv
edges = [(row["node1"], row["node2"]) for _, row in edges_df.iterrows()]

# 创建完整的Graph对象
G = nx.Graph()
G.add_nodes_from(degree_dict.keys())
G.add_edges_from(edges)

# 选择前 300 个节点
top_274_nodes = [
    node
    for node, degree in sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[
        :274
    ]
]

# 创建包含前 300 个节点的子图
subgraph = G.subgraph(top_274_nodes)

# # 获取节点按原始 Degree 值排序
sorted_nodes = sorted(
    {node: degree_dict[node] for node in top_274_nodes}.items(),
    key=lambda x: x[1],
    reverse=True,
)

# # 定义同心圆的层数和每层的节点数
layers = [1, 6, 13, 20, 30, 40, 44, 55, 65]  # 每层节点的数量
radii = [0, 0.8, 1.5, 2.15, 2.8, 3.4, 3.9, 4.4, 4.8]  # 每层圆的半径

# 生成每层圆的节点位置
pos = {}
layer_start = 0
for layer, radius in zip(layers, radii):
    layer_nodes = sorted_nodes[layer_start : layer_start + layer]
    theta = np.linspace(0, 2 * np.pi, len(layer_nodes), endpoint=False)
    for (node, degree), angle in zip(layer_nodes, theta):
        pos[node] = (radius * np.cos(angle), radius * np.sin(angle))
    layer_start += layer

# 获取所有节点的原始 Degree 值
degrees = [degree_dict[node] for node in subgraph.nodes()]

# 生成颜色映射
vmin = np.percentile(degrees, 0.005)
vmax = np.percentile(degrees, 100)
norm = plt.Normalize(vmin, vmax)
cmap = cm.plasma  # 选择一个颜色映射（可以是其他的，如plasma, inferno, magma等）
node_colors = [cmap(norm(degree_dict[node])) for node in subgraph.nodes()]

# 生成节点大小
min_size = 100  # 最小节点大小
max_size = 1000  # 最大节点大小
node_sizes = [
    min_size
    + (max_size - min_size)
    * (degree_dict[node] - min(degrees))
    / (max(degrees) - min(degrees))
    for node in subgraph.nodes()
]

# 生成字体大小
min_font_size = 4  # 最小字体大小
max_font_size = 10  # 最大字体大小
font_sizes = {
    node: min_font_size
    + (max_font_size - min_font_size)
    * (degree_dict[node] - min(degrees))
    / (max(degrees) - min(degrees))
    for node in subgraph.nodes()
}

# 调节线条宽度和透明度,调了，没什么屁用
# edge_width = 0.0001 # 线条宽度
# edge_alpha = 0.0001 # 线条透明度

# 绘制同心圆分布的网络拓扑图
plt.figure(figsize=(6, 6))
nx.draw(
    subgraph,
    pos,
    with_labels=False,
    node_size=node_sizes,
    node_color=node_colors,
    edge_color="gray",
)

# 绘制边
nx.draw_networkx_edges(subgraph, pos, edge_color="gray")

# # 需要调节边的宽度和透明度时用这段绘制，现在还没调成功
# nx.draw_networkx_edges(subgraph, pos, edge_color='gray',width=edge_width,alpha=edge_alpha)

# 手动绘制节点标签并设置字体大小
for node, (x, y) in pos.items():
    plt.text(x, y, s=node, fontsize=font_sizes[node], ha="center", va="center")

plt.title("PPI network top300")
plt.savefig("results/PPI network top300.pdf")
plt.savefig("results/PPI network top300.png")
plt.show()
