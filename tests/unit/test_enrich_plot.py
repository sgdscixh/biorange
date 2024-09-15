import pytest
import pandas as pd
import numpy as np
from plotnine import ggplot
from biorange.enrich_analysis import barplot_go, barplot_kegg


# 创建测试数据
@pytest.fixture
def go_test_data():
    return pd.read_csv("tests/test_data/go_results.csv")


@pytest.fixture
def kegg_test_data():
    return pd.read_csv("tests/test_data/kegg_results.csv")


def test_barplot_go_returns_ggplot(go_test_data):
    plot = barplot_go(go_test_data)
    assert isinstance(plot, ggplot)


def test_barplot_go_correct_number_of_terms(go_test_data):
    plot = barplot_go(go_test_data, shownNumber=2)
    plot_data = plot.data
    assert len(plot_data) == 6  # 6 terms for each of BP, CC, MF


def test_barplot_kegg_returns_ggplot(kegg_test_data):
    plot = barplot_kegg(kegg_test_data)
    assert isinstance(plot, ggplot)


def test_barplot_kegg_correct_number_of_terms(kegg_test_data):
    plot = barplot_kegg(kegg_test_data, shownNumber=7)
    plot_data = plot.data
    assert len(plot_data) == 7


# 测试异常情况
def test_barplot_go_empty_dataframe():
    with pytest.raises(Exception):  # 您可能需要指定更具体的异常类型
        barplot_go(pd.DataFrame())


def test_barplot_kegg_missing_columns(kegg_test_data):
    invalid_data = kegg_test_data.drop(columns=["Overlap"])
    with pytest.raises(Exception):  # 您可能需要指定更具体的异常类型
        barplot_kegg(invalid_data)
