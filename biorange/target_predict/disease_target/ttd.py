import pandas as pd
from biorange.utils.package_fileload import get_data_file_path


class TTDDiseaseScraper:

    def __init__(self, file_path="TTD_combinez_data.csv"):
        self.df = pd.read_csv(get_data_file_path(file_path))

    def search(self, diseases):
        # 如果输入是字符串，则将其转换为包含一个元素的列表
        if isinstance(diseases, str):
            diseases = [diseases]

        # 创建一个空的DataFrame来存储匹配结果
        filtered_df = pd.DataFrame()

        # 遍历输入的diseases列表，进行部分匹配
        for disease in diseases:
            # 使用str.contains进行部分匹配，忽略大小写
            matches = self.df[
                self.df["Disease Entry"].str.contains(disease, case=False, na=False)
            ]
            # 将匹配结果追加到filtered_df中
            filtered_df = pd.concat([filtered_df, matches])

        # 创建一个新的表格，并在第一列增加“data_source”列，内容为“TTD”
        if filtered_df.empty:
            return pd.DataFrame(columns=["disease", "gene_name", "source"])

        new_df = pd.DataFrame()
        new_df["disease"] = filtered_df["Disease Entry"].values
        new_df["gene_name"] = filtered_df["GENENAME"].values
        new_df["source"] = ["TTD"] * len(filtered_df)

        # 将所有非字符串类型的值转换为空字符串
        new_df["gene_name"] = new_df["gene_name"].astype(str)

        # 拆分“gene_name”列中包含多个基因的行
        split_rows = []
        for _, row in new_df.iterrows():
            genes = row["gene_name"].split(";")
            for gene in genes:
                split_rows.append(
                    {
                        "disease": row["disease"],
                        "gene_name": gene.strip(),
                        "source": row["source"],
                    }
                )

        # 创建一个新的DataFrame来存储拆分后的结果
        split_df = pd.DataFrame(split_rows)

        return split_df


ttd_disease_target = TTDDiseaseScraper().search

if __name__ == "__main__":
    searcher = TTDDiseaseScraper()
    result_df = searcher.search("Osteoarthritis")
    result_df.to_csv("results/output2/diseltarge/ttd_output_target.csv", index=False)
    print(result_df)
