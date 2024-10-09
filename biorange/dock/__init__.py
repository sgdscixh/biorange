from biorange.dock.dock_run_diffdock import DiffDockRunner, diffdock_run
from biorange.dock.dock_run_autodock import AutoDockRunner, autodock_run

from biorange.dock.dock_plot import visualize_molecule_complex as plot3d

from biorange.dock.dock_score import smina_score

if __name__ == "__main__":
    x = smina_score(
        ligand_file="results/complex_0/rank1_confidence-0.77.sdf",
        receptor_file="results/6w70.pdb",
    )
    print(x)
