import os
import subprocess

import yaml
from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.tools import ToolClient


class AutoDockRunner:

    def __init__(
        self,
        galaxy_url="https://usegalaxy.org/",
        api_key="671a2f359458e80d29a18ffa0a8683f1",  # TODO:全局配置
    ):
        self.galaxy_url = galaxy_url
        self.api_key = api_key
        self.gi = GalaxyInstance(galaxy_url, api_key)
        self.tool_client = ToolClient(self.gi)

    def create_history(self, name):
        history = self.gi.histories.create_history(name=name)
        return history["id"]

    def upload_file(self, file_path, file_type, history_id):
        upload_result = self.tool_client.upload_file(
            file_path, history_id, file_type=file_type
        )
        return upload_result["outputs"][0]["id"]

    def run_planemo(
        self,
        outdir,
        pdb_galaxy_id,
        sdf_galaxy_id,
        history_name="AutoDock Test Run",
        tags="autodock,test-run",
        tool_id="aee83f8c6aa459ef",
        planemo_executable="planemo",
    ):
        job_data = {
            "pdb_file": {"class": "File", "galaxy_id": pdb_galaxy_id},
            "sdf_file": {"class": "File", "galaxy_id": sdf_galaxy_id},
        }

        job_file_path = os.path.join(outdir, "job.yml")
        os.makedirs(outdir, exist_ok=True)
        with open(job_file_path, "w", encoding="utf-8") as job_file:
            yaml.dump(job_data, job_file, default_flow_style=False)

        command = [
            planemo_executable,
            "run",
            tool_id,
            job_file_path,
            "--engine",
            "external_galaxy",
            "--galaxy_url",
            self.galaxy_url,
            "--galaxy_user_key",
            self.api_key,
            "--download_outputs",
            "--outdir",
            outdir,
            "--history_name",
            history_name,
            "--tags",
            tags,
        ]

        try:
            subprocess.run(command, check=True)
            print("Planemo run completed successfully.")
            os.remove(job_file_path)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred during planemo run: {e}")

    def runflow(self, pdb_file_path, sdf_file_path, outdir="."):
        # 创建新的历史记录
        history_id = self.create_history("测试全流程")

        # 上传文件并获取ID
        pdb_file_id = self.upload_file(pdb_file_path, "pdb", history_id)
        sdf_file_id = self.upload_file(sdf_file_path, "sdf", history_id)

        # 运行Planemo
        self.run_planemo(
            outdir=outdir,
            pdb_galaxy_id=pdb_file_id,
            sdf_galaxy_id=sdf_file_id,
            history_name="autodock",
            tags="autodock,biorange",
        )

        print("Workflow completed successfully.")
        # 清理文件
        try:
            os.remove("tool_test_output.html")
            os.remove("tool_test_output.json")
        except FileNotFoundError:
            pass


autodock_run = AutoDockRunner().runflow

# 使用示例
if __name__ == "__main__":
    workflow = AutoDockRunner()
    pdb_file = "/home/liuyan/projects/package/biorange/biorange/dock/todo/6w70.pdb"
    sdf_file = (
        "/home/liuyan/projects/package/biorange/biorange/dock/todo/6w70_ligand.sdf"
    )
    outdir = "/home/liuyan/projects/package/biorange/biorange/dock/todo/"

    workflow.runflow(pdb_file, sdf_file, outdir)
