"""
模块: __init__.py
描述: 'biorange' 包的初始化模块。
"""

from pyfiglet import figlet_format
from termcolor import colored
from . import component, enrich_analysis, ppi, target_predict, utils, venn


# 尝试获取已安装包的版本号
try:
    from importlib.metadata import version as get_version

    __version__ = get_version("biorange")
except ImportError:
    __version__ = "0.0.0"  # 如果获取版本号失败，则使用默认版本号

from pyfiglet import figlet_format


# 生成斜体的 "BIORANGE" ASCII 艺术字
ascii_art = figlet_format("BioRange", font="slant")
colored_art = colored(ascii_art, color="cyan", attrs=["bold"])

# 打印带有边框的艺术字和版本号
print(colored_art.center(50))
print(colored(f"Version: {__version__}".center(50), color="yellow", attrs=["bold"]))

# 控制导出的东西，隐藏细节
## 只有被all纳入的才会被导出，不设置__all__的话默认暴露所有不以下划线开头的变量和函数
__all__ = [
    "component",
    "enrich_analysis",
    "ppi",
    "target_predict",
    "utils",
    "venn",
]


def init_print():
    """初始化日志记录器。"""
    pass
