import pytest
import pandas as pd
from typing import List
from unittest.mock import patch, MagicMock
from biorange.enrich_analysis import (
    enrich_gokegg,
)


# 很明显这里包装了一个输入数据用来测试
@pytest.fixture
def sample_gene_list():
    return [
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


@pytest.fixture
def mock_enrichr_results_all():
    mock_results = pd.read_csv("tests/test_data/all_results.csv")
    return mock_results


@pytest.fixture
def mock_enrichr_results_go():
    mock_results = pd.read_csv("tests/test_data/go_results.csv")
    return mock_results


@pytest.fixture
def mock_enrichr_results_kegg():
    mock_results = pd.read_csv("tests/test_data/kegg_results.csv")
    return mock_results


@patch("gseapy.enrichr")
def test_enrich_gokegg_all(mock_enrichr, sample_gene_list, mock_enrichr_results_all):
    mock_enrichr.return_value = MagicMock(results=mock_enrichr_results_all)

    result = enrich_gokegg(sample_gene_list, analysis_type="all")

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 666
    assert all(result["Adjusted P-value"] < 0.05)
    assert mock_enrichr.call_count == 2


@patch("gseapy.enrichr")
def test_enrich_gokegg_go(mock_enrichr, sample_gene_list, mock_enrichr_results_go):
    mock_enrichr.return_value = MagicMock(
        results=mock_enrichr_results_go[
            mock_enrichr_results_go["Gene_set"] == "GO_Biological_Process_2021"
        ]
    )

    result = enrich_gokegg(sample_gene_list, analysis_type="go")

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 274
    assert all(result["Gene_set"] == "GO_Biological_Process_2021")
    assert mock_enrichr.call_count == 1


@patch("gseapy.enrichr")
def test_enrich_gokegg_kegg(mock_enrichr, sample_gene_list, mock_enrichr_results_kegg):
    mock_enrichr.return_value = MagicMock(
        results=mock_enrichr_results_kegg[
            mock_enrichr_results_kegg["Gene_set"] == "KEGG_2021_Human"
        ]
    )

    result = enrich_gokegg(sample_gene_list, analysis_type="kegg")

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 25
    assert all(result["Gene_set"] == "KEGG_2021_Human")
    assert mock_enrichr.call_count == 1
