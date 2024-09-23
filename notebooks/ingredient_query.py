import pandas as pd
from biorange.component import TCMSPComponentLocalScraper

# 初始化搜索客户端
search_client = TCMSPComponentLocalScraper(use_remote=True)

# 草药列表
herb_list = ["半夏", "人参", "黄芩", "干姜", "甘草", "大枣", "黄连"]


def fetch_herb_components(herb_list):
    dataframes = []
    for herb in herb_list:
        df = search_client.search_herb(herb)  # 假设返回的是一个 DataFrame
        df["herb"] = herb  # 添加一列来标识草药名称
        dataframes.append(df)
    return dataframes


# 获取查询结果
dataframes = fetch_herb_components(herb_list)

# 合并所有 DataFrame
combined_df = pd.concat(dataframes, ignore_index=True)

# 保存到 CSV 文件
combined_df.to_csv("herb_components.csv", index=False)

print("数据已保存到 herb_components.csv")
