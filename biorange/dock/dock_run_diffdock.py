import os
import subprocess
from pathlib import Path
from biorange.utils.package_fileload import get_data_file_path


from biorange.logger import get_logger

logger = get_logger(__name__)


class DiffDockRunner:
    def __init__(self, docker_image, entrypoint="/bin/bash"):
        self.docker_image = docker_image
        self.entrypoint = entrypoint

    def run_inference(
        self,
        protein_path,
        ligand_path,
        samples_per_complex=10,
        output_dir="results",
        config_path=get_data_file_path("config-dock.yaml"),
    ):
        """
        Run DiffDock inference using Docker.

        Args:
            protein_path (str): Path to the protein file.
            ligand_path (str): Path to the ligand file.
            config_path (str): Path to the configuration file.
            samples_per_complex (int): Number of samples per complex.
            output_dir (str): Directory to store the output.

        Raises:
            FileNotFoundError: If any input file is not found.
            subprocess.CalledProcessError: If the Docker command fails.
        """
        self._validate_inputs(
            protein_path, ligand_path, config_path, samples_per_complex, output_dir
        )

        docker_cmd = self._build_docker_command(
            protein_path, ligand_path, config_path, samples_per_complex, output_dir
        )
        # 修改容器外部的权限
        logger.info("本分析在容器中进行，重新运行需要清空文件夹")
        os.makedirs(f"{output_dir}/complex_0")
        subprocess.run(["chmod", "-R", "777", output_dir], check=True)

        try:
            logger.info("Starting DiffDock inference...")
            result = subprocess.run(
                docker_cmd, check=True, text=True, capture_output=True
            )
            logger.info("DiffDock inference completed successfully.")
            logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            logger.error(f"Docker command failed: {e}")
            logger.error(f"Error output: {e.stderr}")
            raise

    def _validate_inputs(
        self, protein_path, ligand_path, config_path, samples_per_complex, output_dir
    ):
        for path in [protein_path, ligand_path, config_path]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"File not found: {path}")

        if not isinstance(samples_per_complex, int) or samples_per_complex <= 0:
            raise ValueError("samples_per_complex must be a positive integer")

        Path(output_dir).mkdir(parents=True, exist_ok=True)

    def _build_docker_command(
        self, protein_path, ligand_path, config_path, samples_per_complex, output_dir
    ):
        # Convert all paths to absolute paths
        protein_path = os.path.abspath(protein_path)
        ligand_path = os.path.abspath(ligand_path)
        config_path = os.path.abspath(config_path)
        output_dir = os.path.abspath(output_dir)

        # Define container paths
        container_protein_dir = "/mnt/protein"
        container_ligand_dir = "/mnt/ligand"
        container_config_dir = "/mnt/config"
        container_output_dir = "/mnt/output"

        # Build volume mappings
        volume_mappings = [
            f"-v{os.path.dirname(protein_path)}:{container_protein_dir}",
            f"-v{os.path.dirname(ligand_path)}:{container_ligand_dir}",
            f"-v{os.path.dirname(config_path)}:{container_config_dir}",
            f"-v{output_dir}:{container_output_dir}",
        ]

        # Build the Docker command
        docker_cmd = [
            "docker",
            "run",
            "--rm",
            *volume_mappings,
            "--entrypoint",
            self.entrypoint,
            self.docker_image,
            "-c",
            (
                f"micromamba run -n diffdock python inference.py "
                f"--protein_path={container_protein_dir}/{os.path.basename(protein_path)} "
                f"--ligand={container_ligand_dir}/{os.path.basename(ligand_path)} "
                f"--config={container_config_dir}/{os.path.basename(config_path)} "
                f"--samples_per_complex={samples_per_complex} "
                f"--out_dir={container_output_dir}"
            ),
        ]

        logger.info(f"Docker command: {' '.join(docker_cmd)}")
        return docker_cmd


diffdock_run = DiffDockRunner(
    docker_image="august777/diffdocker:full-0910", entrypoint="/bin/bash"
).run_inference

if __name__ == "__main__":
    # 设置参数
    protein_path = "results/6w70.pdb"
    ligand_path = "results/6w70_ligand.sdf"
    docker_image = "august777/diffdocker:full-0910"

    runner = DiffDockRunner(docker_image=docker_image)

    runner.run_inference(
        protein_path=protein_path,
        ligand_path=ligand_path,
    )
