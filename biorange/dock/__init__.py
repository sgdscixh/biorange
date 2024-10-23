"""分子对接模块"""

# 预处理
from biorange.dock.dock_clean import fix_pdb as dock_fix_protein

# 分子对接软件封装
from biorange.dock.dock_run_diffdock import DiffDockRunner, diffdock_run
from biorange.dock.dock_run_autodock import AutoDockRunner, autodock_run

# 绘图
from biorange.dock.dock_plot import visualize_molecule_complex as plot3d
from biorange.dock.dock_plot_interaction import (
    ProteinLigandAnalyzer,
    dock_plot_2d,
    dock_plot_3d,
    dock_interaction_df,
)

# 评分
from biorange.dock.dock_score import smina_score

if __name__ == "__main__":
    x = smina_score(
        ligand_file="results/complex_0/rank1_confidence-0.77.sdf",
        receptor_file="results/6w70.pdb",
    )
    print(x)
