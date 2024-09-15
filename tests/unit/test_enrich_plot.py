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


# 测试就是一种限制~~限制后面的人别瞎改，也是一种期待，期待代码具有怎样的表现
# 然后写出符合这个期待的代码。
# 因此存在两种模式  一种是先写代码然后写测试，一种是先写期待的测试结果，再写代码去通过测试


## 这里我期待透明度设置的标签是log10(P-value)


## 成功啦
## 接下来我开始模拟一不小心修改了代码 "Alpha scale label should be '-log10(P-value)'"
# E       AssertionError: Alpha scale label should be '-log10(P-value)'
# E       assert '-log10(P-哈哈value)' == '-log10(P-value)'


## 果然报错啦，然后我就收到了提醒。必须修改这个代码 通过测试


## 如何如何~~测试代码主要是模拟用户使用的情况嘛
## 对的这是他一半的目的 和，一半是模拟用户使用，术语叫做功能测试functional test  另一部分测试代码有无报错异常学名单元测试 两着难以区分oo

## 所以总结完成测试流程就是
# 心中有期待，然后告诉ai，然后生成代码，和测试代码 哈哈 悟了吗？okk,开始我想着直接在源代码上直接改也可，但是这种包是汇聚了很多功能的，用户不是直接在源代码上用
# 啊对对对，源代码上测试也很流行，但是不够~~这样是最佳实践！
# 所以看不懂测试代码不重要  都是ai写的 心中期待什么功能很重要
# 比如你不期待重复值，这些都应该通过测试题限制  而不是手动搓搓~~
# 测试代码不重要，主要是是能实现功能对的~证明源代码是可以实现那些功能的完全正确！！举一反三~~b(￣▽￣)d舒服了吗okk你下班哦哦哦
# 还有就是人类眼睛注意不到那么多情况~~不小心改错代码经常有~~测试可以固定行为 每次提交的时候看看测试还是不是正常的


# 另外还有测试覆盖率的概念，就是边边角角都要测试到 那么系统就是无敌的了~~文件的一批 溜啦溜啦哦空空ok 你继续玩玩~~
def test_go_bar_plot_alpha_label():
    from plotnine.scales import scale_alpha_continuous

    # 创建测试数据
    df = pd.DataFrame(
        {
            "Gene_set": ["BP", "CC", "MF"] * 5,
            "Description": [f"Term_{i}" for i in range(15)],
            "Count": np.random.randint(10, 100, 15),
            "P-value": np.random.uniform(0.001, 0.05, 15),
        }
    )

    # 生成图形
    p = barplot_go(df)

    # 直接访问图形对象中的透明度刻度标签
    alpha_scale = next(
        (scale for scale in p.scales if isinstance(scale, scale_alpha_continuous)), None
    )

    # 检查透明度刻度标签是否设置为 "-log10(P-value)"
    assert alpha_scale is not None, "Alpha scale should be present in the plot."
    assert (
        alpha_scale.name == "-log10(P-value)"
    ), "Alpha scale label should be '-log10(P-value)'"
