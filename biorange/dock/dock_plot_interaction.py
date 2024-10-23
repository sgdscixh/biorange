import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
from typing import Union, List, Optional, Tuple
import pandas as pd


class ProteinLigandAnalyzer:
    """蛋白质-配体相互作用分析器

    用于分析蛋白质与配体之间的相互作用，支持2D和3D可视化。
    每个方法都可以独立运行，无需考虑调用顺序。
    """

    def __init__(self):
        self.protein_mol = None
        self.ligand_mol = None
        self.fingerprint = None

    def _ensure_structures(
        self,
        protein_file: Optional[str] = None,
        ligand_file: Optional[str] = None,
        protein_method: str = "mda",
    ) -> Tuple[plf.Molecule, plf.Molecule]:
        """确保结构已加载，如果没有则加载

        Args:
            protein_file: 可选的蛋白质文件路径
            ligand_file: 可选的配体文件路径
            protein_method: 蛋白质加载方法

        Returns:
            protein_mol, ligand_mol: 处理后的分子对象
        """
        protein_mol = self.protein_mol
        ligand_mol = self.ligand_mol

        if protein_file is not None or protein_mol is None:
            if protein_file is None:
                raise ValueError("Protein file is required")
            if protein_method == "mda":
                u = mda.Universe(protein_file)
                u.atoms.guess_bonds(vdwradii={"H": 1.05, "O": 1.48})
                protein_mol = plf.Molecule.from_mda(u)
            elif protein_method == "rdkit":
                rdkit_prot = Chem.MolFromPDBFile(protein_file, removeHs=False)
                protein_mol = plf.Molecule(rdkit_prot)
            self.protein_mol = protein_mol

        if ligand_file is not None or ligand_mol is None:
            if ligand_file is None:
                raise ValueError("Ligand file is required")
            ligand_mol = plf.sdf_supplier(ligand_file)[0]
            self.ligand_mol = ligand_mol

        return protein_mol, ligand_mol

    def _ensure_fingerprint(
        self, protein_mol: plf.Molecule, ligand_mol: plf.Molecule, count: bool = True
    ) -> plf.Fingerprint:
        """确保指纹已生成，如果没有则生成

        Args:
            protein_mol: 蛋白质分子对象
            ligand_mol: 配体分子对象
            count: 是否计数所有可能的相互作用组合

        Returns:
            fingerprint: 生成的指纹对象
        """
        fingerprint = plf.Fingerprint(count=count)
        fingerprint.run_from_iterable([ligand_mol], protein_mol)
        self.fingerprint = fingerprint
        return fingerprint

    def get_interaction_df(
        self,
        protein_file: Optional[str] = None,
        ligand_file: Optional[str] = None,
        protein_method: str = "mda",
        count: bool = True,
    ) -> pd.DataFrame:
        """获取相互作用数据框

        Args:
            protein_file: 可选的蛋白质文件路径
            ligand_file: 可选的配体文件路径
            protein_method: 蛋白质加载方法
            count: 是否计数所有可能的相互作用组合

        Returns:
            包含相互作用信息的DataFrame
        """
        protein_mol, ligand_mol = self._ensure_structures(
            protein_file, ligand_file, protein_method
        )
        fingerprint = self._ensure_fingerprint(protein_mol, ligand_mol, count)
        return fingerprint.to_dataframe()

    def visualize_2d(
        self,
        protein_file: Optional[str] = None,
        ligand_file: Optional[str] = None,
        protein_method: str = "mda",
        display_all: bool = True,
        count: bool = True,
        save: Optional[str] = None,
    ):
        """生成2D交互网络图

        Args:
            protein_file: 可选的蛋白质文件路径
            ligand_file: 可选的配体文件路径
            protein_method: 蛋白质加载方法
            display_all: 是否显示所有可能的相互作用
            count: 是否计数所有可能的相互作用组合
            save: 可选的保存路径

        Returns:
            2D网络图视图对象
        """
        protein_mol, ligand_mol = self._ensure_structures(
            protein_file, ligand_file, protein_method
        )
        fingerprint = self._ensure_fingerprint(protein_mol, ligand_mol, count)
        view_2d = fingerprint.plot_lignetwork(
            ligand_mol, kind="frame", frame=0, display_all=display_all
        )

        if save:
            with open(save, "w", encoding="utf-8") as f:
                f.write(view_2d.data)

        return view_2d

    def visualize_3d(
        self,
        protein_file: Optional[str] = None,
        ligand_file: Optional[str] = None,
        protein_method: str = "mda",
        display_all: bool = False,
        count: bool = True,
        save: Optional[str] = None,
    ):
        """生成3D结构视图

        Args:
            protein_file: 可选的蛋白质文件路径
            ligand_file: 可选的配体文件路径
            protein_method: 蛋白质加载方法
            display_all: 是否显示所有可能的相互作用
            count: 是否计数所有可能的相互作用组合
            save: 可选的保存路径

        Returns:
            3D结构视图对象
        """
        protein_mol, ligand_mol = self._ensure_structures(
            protein_file, ligand_file, protein_method
        )
        fingerprint = self._ensure_fingerprint(protein_mol, ligand_mol, count)
        view_3d = fingerprint.plot_3d(
            ligand_mol, protein_mol, frame=0, display_all=display_all
        )

        if save:
            view_3d.write_html(save)

        return view_3d


dock_plot_2d = ProteinLigandAnalyzer().visualize_2d
dock_plot_3d = ProteinLigandAnalyzer().visualize_3d
dock_interaction_df = ProteinLigandAnalyzer().get_interaction_df

if __name__ == "__main__":
    protein_file = str(
        plf.datafiles.datapath / "vina" / "rec.pdb"
    )  # plf.datafiles.datapath / "vina" / "rec.pdb"
    ligand_file = str(
        plf.datafiles.datapath / "vina" / "vina_output.sdf"
    )  # plf.datafiles.datapath / "vina" / "vina_output.sdf"
    # 生成并保存2D图
    view_2d = ProteinLigandAnalyzer().visualize_2d(
        protein_file, ligand_file, save="results/dock/saving.html"
    )
    # 生成并保存3D图
    view_3d = ProteinLigandAnalyzer().visualize_3d(
        protein_file, ligand_file, save="results/dock/dd.html"
    )

# 这是简化版的绘图
# 真正的英文原版在https://github.com/chemosim-lab/ProLIF/blob/master/docs/notebooks/pdb.ipynb
# 这个仓库很不错，里面还有动力学和对接的专业教程
# 真的成功了嘿嘿
