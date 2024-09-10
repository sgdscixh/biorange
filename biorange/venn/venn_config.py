# Define configuration using a dataclass for clarity and type safety
from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class VennPlotConfig:
    font_family: str = "Times New Roman"
    figsize: Tuple[int, int] = (5, 3)
    dpi: int = 110
    set_colors: List[str] = field(
        default_factory=lambda: ["#0073C2FF", "#EFC000FF", "#CD534CFF"]
    )
    alpha: float = 0.8
    linestyle: str = "--"
    linewidth: float = 1.7
    circle_color: str = "black"
    title_fontsize: int = 12
    label_fontsize: int = 10
    label_fontweight: str = "normal"
    title_fontweight: str = "bold"
    title_color: str = "black"
    title_style: str = "normal"

    def update(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise KeyError(f"Invalid configuration key: {key}")
