# 分割

from openbabel import pybel
from pathlib import Path


# 将 SDF 文件拆分为单独的分子文件
def split_sdf_file(sdf_path):
    """
    Split an SDF file into seperate files for each molecule.
    Each file is named with consecutive numbers.

    Parameters
    ----------
    sdf_path: str or pathlib.Path
        Path to SDF file that should be split.
    """
    sdf_path = Path(sdf_path)
    stem = sdf_path.stem
    parent = sdf_path.parent
    molecules = pybel.readfile("sdf", str(sdf_path))
    # 将每个分子写入单独的 SDF 文件
    for i, molecule in enumerate(molecules, 1):
        molecule.write("sdf", str(parent / f"{stem}_{i}.sdf"), overwrite=True)
    return


# 使用示例
split_sdf_file(
    "/home/liuyan/projects/package/biorange/biorange/dock/todo/ligand1_docked__3e8dbfee-2304-4636-a3ec-8b8e57975564"
)
