import gseapy as gp
import pandas as pd
from typing import List, Literal


def perform_enrichment_analysis(
    gene_list: List[str],
    analysis_type: Literal["go", "kegg", "all"] = "all",
    organism: str = "human",
    kegg_database: str = "KEGG_2021_Human",
    go_databases: List[str] = None,
    p_value_cutoff: float = 0.05,
) -> pd.DataFrame:
    """
    对给定的基因列表执行GO和/或KEGG富集分析。

    参数:
        gene_list (List[str]): 要分析的基因符号列表。
        analysis_type (Literal['go', 'kegg', 'all']): 要执行的分析类型。默认为'all'。
        organism (str): 分析所针对的物种。默认为 "human"。
        kegg_database (str): 使用的KEGG数据库。默认为 "KEGG_2021_Human"。
        go_databases (List[str]): 使用的GO数据库列表。如果为None，则使用默认的GO数据库。
        p_value_cutoff (float): 显著性的p值阈值。默认为0.05。

    返回:
        pd.DataFrame: 包含富集分析结果的DataFrame。
    """
    # 如果未指定GO数据库，则使用默认的GO数据库
    if go_databases is None:
        go_databases = [
            "GO_Biological_Process_2021",
            "GO_Molecular_Function_2021",
            "GO_Cellular_Component_2021",
        ]

    results = pd.DataFrame()  # 初始化空的DataFrame来存储结果

    # 如果分析类型是GO或全部
    if analysis_type in ["go", "all"]:
        # 使用gseapy进行GO富集分析
        enr_go = gp.enrichr(
            gene_list=gene_list,
            gene_sets=go_databases,
            organism=organism,
            cutoff=1,  # 设置为1以获取所有结果，稍后会过滤
            no_plot=True,
            outdir=None,
        )
        go_results = enr_go.results
        results = pd.concat([results, go_results])  # 将GO结果添加到总结果中

    # 如果分析类型是KEGG或全部
    if analysis_type in ["kegg", "all"]:
        # 使用gseapy进行KEGG富集分析
        enr_kegg = gp.enrichr(
            gene_list=gene_list,
            gene_sets=[kegg_database],
            organism=organism,
            cutoff=1,  # 设置为1以获取所有结果，稍后会过滤
            no_plot=True,
            outdir=None,
        )
        kegg_results = enr_kegg.results
        results = pd.concat([results, kegg_results])  # 将KEGG结果添加到总结果中

    # 根据p值阈值过滤结果
    results = results[results["Adjusted P-value"] < p_value_cutoff]

    # 根据调整后的p值对结果排序
    results = results.sort_values("Adjusted P-value")

    return results


# 示例用法
if __name__ == "__main__":
    gene_list = [
        "MLH1",
        "ECM10",
        "RLI1",
        "SSB1",
        "SSB2",
        "MSH2",
        "STAT3",
        "IL6",
        "S100A9",
        "S100A8",
        "ARG1",
        "NOS2",
    ]

    # 执行GO和KEGG分析
    all_results = perform_enrichment_analysis(gene_list, analysis_type="all")
    print("GOKEGG集分析结果:")
    print(all_results.head(10))

    # # 仅执行GO分析
    # go_results = perform_enrichment_analysis(gene_list, analysis_type="go")
    # print("\nGO富集分析结果:")
    # print(go_results.head(10))

    # # 仅执行KEGG分析 默认就是分析go+kegg
    # kegg_results = perform_enrichment_analysis(gene_list, analysis_type="kegg")
    # print("\nKEGG富集分析结果:")
    # print(kegg_results.head(10))
