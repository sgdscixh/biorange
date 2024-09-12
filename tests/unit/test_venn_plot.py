import pytest
from unittest.mock import patch, MagicMock
import matplotlib.pyplot as plt
from biorange.venn.venn_plot import VennPlotter, VennPlotConfig

# Test data
TCMSP = ["A", "B", "C", "D"]
Chembl = ["B", "C", "E", "F"]
STITCH = ["C", "D", "E", "G"]
TTD = ["H", "I", "J", "D", "E", "F", "G"]


@pytest.fixture
def venn_plotter():
    return VennPlotter()


def test_venn_plotter_initialization():
    plotter = VennPlotter()
    assert isinstance(plotter.config, VennPlotConfig)
    assert plt.rcParams["font.family"] == [plotter.config.font_family]


def test_venn_plotter_custom_config():
    custom_config = VennPlotConfig(figsize=(8, 6), dpi=150, title_fontsize=16)
    plotter = VennPlotter(config=custom_config)
    assert plotter.config.figsize == (8, 6)
    assert plotter.config.dpi == 150
    assert plotter.config.title_fontsize == 16


@patch("pandas.DataFrame.to_csv")
@patch("matplotlib.pyplot.savefig")
@patch("matplotlib.pyplot.show")
def test_plot_venn_2_groups(mock_to_csv, mock_show, mock_savefig, venn_plotter):
    venn_plotter.plot_venn(
        [TCMSP, Chembl], ["TCMSP", "Chembl"], "Test Venn Diagram", "test_venn_2.png"
    )
    mock_savefig.assert_called_once()
    mock_show.assert_called_once()


@patch("matplotlib.pyplot.savefig")
@patch("matplotlib.pyplot.show")
@patch("pandas.DataFrame.to_csv")
def test_plot_venn_3_groups(mock_to_csv, mock_show, mock_savefig, venn_plotter):
    venn_plotter.plot_venn(
        [TCMSP, Chembl, STITCH],
        ["TCMSP", "Chembl", "STITCH"],
        "Test Venn Diagram",
        "test_venn_3.png",
    )
    mock_savefig.assert_called_once()
    mock_to_csv.assert_called_once()
    mock_show.assert_called_once()


def test_plot_venn_invalid_groups(venn_plotter):
    with pytest.raises(ValueError, match="Only 2 or 3 groups are supported."):
        venn_plotter.plot_venn(
            [TCMSP, Chembl, STITCH, ["H", "I"]],
            ["TCMSP", "Chembl", "STITCH", "Extra"],
            "Invalid Venn Diagram",
        )


@patch("pandas.DataFrame.to_csv")
@patch("matplotlib.pyplot.savefig")
@patch("matplotlib.pyplot.show")
def test_plot_venn_dynamic_config_update(
    mock_to_csv, mock_show, mock_savefig, venn_plotter
):
    venn_plotter.plot_venn(
        [TCMSP, Chembl],
        ["TCMSP", "Chembl"],
        "Updated Venn Diagram",
        "test_updated_venn.png",
        title_fontsize=18,
        label_fontsize=14,
    )
    assert venn_plotter.config.title_fontsize == 18
    assert venn_plotter.config.label_fontsize == 14
    mock_savefig.assert_called_once()
    mock_show.assert_called_once()


##############################################################################
# 上面的不用看，下面的接近使用逻辑
@patch("os.makedirs")
def test_plot_venn_create_output_directory_2_groups(mock_makedirs, venn_plotter):
    venn_plotter.plot_venn(
        [TCMSP, Chembl],
        ["TCMSP", "Chembl"],
        "Test 2 Groups",
        "test_create_dir_2_groups.png",
    )
    mock_makedirs.assert_called_once_with("results/venn", exist_ok=True)


@patch("os.makedirs")
def test_plot_venn_create_output_directory_3_groups(mock_makedirs, venn_plotter):
    venn_plotter.plot_venn(
        [TCMSP, Chembl, STITCH],
        ["TCMSP", "Chembl", "STITCH"],
        "Test 3 Groups",
        "test_create_dir_3_groups.png",
    )
    mock_makedirs.assert_called_once_with("results/venn", exist_ok=True)


@patch("os.makedirs")
def test_plot_venn_create_output_directory_33_groups(mock_makedirs, venn_plotter):
    venn_plotter.plot_venn(
        groups=[TCMSP, Chembl, TTD],  # TTD 没给实例列表 你可以新增
        labels=["TCMSP", "Chembl", "TTD"],
        title="Test 33 Groups",  # pytest 会自动执行test文件夹下所有test_开头的函数。类
        filename="test_create_dir_33_groups.png",
    )
    mock_makedirs.assert_called_once_with("results/venn", exist_ok=True)
