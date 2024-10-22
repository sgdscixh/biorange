import re
import subprocess


def smina_score(receptor_file, ligand_file):
    # 构建 smina 命令
    command = ["smina", "-r", receptor_file, "-l", ligand_file, "--score_only"]

    try:
        # 运行 smina 命令
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        output = result.stdout

        # 提取 Affinity 和 Intramolecular energy
        affinity_match = re.search(r"Affinity:\s*([-\d.]+)", output)
        intra_energy_match = re.search(r"Intramolecular energy:\s*([-\d.]+)", output)

        affinity = float(affinity_match.group(1)) if affinity_match else None
        intra_energy = (
            float(intra_energy_match.group(1)) if intra_energy_match else None
        )

        return {
            "affinity": affinity,
            "intramolecular_energy": intra_energy,
            "all": output,
        }

    except subprocess.CalledProcessError as e:
        print(f"Error running smina: {e}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


if __name__ == "__main__":

    # 使用示例
    # 先运行dock_run_diffdock.py
    receptor_file = "/home/liuyan/projects/package/biorange/biorange/dock/todo/AF-P31749-F1-model_v4.pdb"
    ligand_file = "/home/liuyan/projects/package/biorange/results/complex_0/rank1_confidence-0.18.sdf"
    # receptor_file = "/home/liuyan/projects/package/biorange/notebooks/fixed_output.pdb"
    # ligand_file = "/home/liuyan/projects/package/biorange/results/complex_0/rank1_confidence-0.18.sdf"

    results = smina_score(receptor_file, ligand_file)
    print(results)
