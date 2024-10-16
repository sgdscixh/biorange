from typing import List, Optional
import pandas as pd
import requests
from biorange.utils.package_fileload import get_data_file_path


class ChemblTargetPredictor:
    """
    用于从ChEMBL获取目标预测的类。2019版离线chembl,可无限使用。
    """

    def __init__(self):
        self.url = " http://127.0.0.1:8080"
        self.headers = {"Content-Type": "application/json"}

    def predict_single(self, smiles: str) -> Optional[pd.DataFrame]:
        """
        为单个SMILES字符串获取目标预测。

        参数:
        smiles (str): 输入的SMILES字符串

        返回:
        Optional[pd.DataFrame]: 包含预测结果的DataFrame,如果失败则返回None
        """
        try:
            response = requests.post(
                self.url, json={"smiles": smiles}, headers=self.headers
            )
            response.raise_for_status()
            predictions = response.json()

            results = []
            for prediction in predictions:
                result = {
                    "smiles": smiles,
                    "targetChemblId": prediction.get("target_chemblid"),
                    "organism": prediction.get("organism"),
                    "prefName": prediction.get("pref_name"),
                    "confidence70": prediction.get("70%"),
                    "confidence80": prediction.get("80%"),
                    "confidence90": prediction.get("90%"),
                    "threshold": prediction.get("threshold"),
                }
                results.append(result)

            return pd.DataFrame(results)

        except requests.exceptions.RequestException as e:
            print(f"请求失败: {e}")
            return None

    def predict_multiple(self, smiles_list: List[str]) -> pd.DataFrame:
        """
        为多个SMILES字符串获取目标预测并合并结果。

        参数:
        smiles_list (List[str]): SMILES字符串列表

        返回:
        pd.DataFrame: 包含所有预测结果的合并DataFrame
        """
        all_predictions = []
        for smiles in smiles_list:
            df = self.predict_single(smiles)
            if df is not None:
                all_predictions.append(df)

        if all_predictions:
            return pd.concat(all_predictions, ignore_index=True)
        else:
            return pd.DataFrame()

    def additional_processing(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        对预测结果进行额外处理，包括与内置数据表匹配。

        参数:
        df (pd.DataFrame): 包含预测结果的DataFrame

        返回:
        pd.DataFrame: 处理后的DataFrame
        """
        # 读取内置数据文件
        try:
            internal_data_path = get_data_file_path(
                "chembl_uniport_gene25.csv"
            )  # 没问题啊但是他说没有这个文件 在哪里没有呀
            internal_data = pd.read_csv(internal_data_path)
        except FileNotFoundError:
            print("内置数据文件未找到")
            return df

        # 合并预测结果与内置数据表
        merged_df = pd.merge(df, internal_data, on="targetChemblId", how="left")
        return merged_df

    def search_smiles(self, smiles: List[str] | str):
        """
        根据SMILES字符串搜索目标预测。
        """
        if isinstance(smiles, str):
            smiles = [smiles]
        temp_df = self.predict_multiple(smiles)
        return self.additional_processing(temp_df)


if __name__ == "__main__":
    df = pd.read_csv(
        "/home/liuyan/projects/package/biorange/biorange/data/sup_smiles.csv"
    )
    smiles = df["smiles"].tolist()
    predictors = ChemblTargetPredictor()
    res = predictors.search_smiles(smiles)
    res.to_csv("./CHEMBL_25_all_sup.csv")
    print(res.columns.tolist())
