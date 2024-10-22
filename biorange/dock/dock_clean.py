# 给从alphafold下载是蛋白质pbd文件加氢
from pdbfixer import PDBFixer
from openmm.app import PDBFile


# openmm的深度学习方法加氢 先进
def fix_pdb(input_path, output_path, ph=7.4):
    """
    修复PDB文件的主要功能

    Parameters:
    -----------
    input_path : str
        输入PDB文件路径
    output_path : str
        输出PDB文件路径
    ph : float
        设置pH值，默认7.0
    """
    # 初始化fixer
    fixer = PDBFixer(input_path)

    # 查找和添加缺失的残基
    fixer.findMissingResidues()

    # 查找缺失的原子
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # 添加缺失的氢原子
    fixer.addMissingHydrogens(ph)

    # 保存修复后的结构
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_path, "w"))

    print(f"Structure has been fixed and saved to {output_path}")


# 使用示例
fix_pdb(
    "/home/liuyan/projects/package/biorange/biorange/dock/todo/AF-P05177-F1-model_v4.pdb",
    "fixed_output2.pdb",
)
