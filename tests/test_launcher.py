import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class OneClickLauncherTests(unittest.TestCase):
    def test_windows_launcher_has_required_one_click_steps(self):
        text = (ROOT / "ONT_Plasmid_Analyzer_One_Click.bat").read_text(
            encoding="utf-8"
        )
        for expected in (
            "182.198.164.21",
            "MCET03",
            "18502",
            "MobaXterm",
            "-newtab",
            "ont_one_click.sh",
            "Invoke-WebRequest",
            "start \"\" \"%UI_URL%\"",
        ):
            self.assertIn(expected, text)
        self.assertNotIn("“", text)
        self.assertNotIn("”", text)

    def test_server_launcher_uses_expected_environment_and_data_root(self):
        text = (ROOT / "ont_one_click.sh").read_text(encoding="utf-8")
        for expected in (
            "/home/mcet/anaconda3",
            '$HOME/conda_envs/NGS_ONT_env',
            'PROJECT_DIR="${ONT_PROJECT_DIR:-$SCRIPT_DIR}"',
            "UI_ADDRESS=\"${ONT_UI_ADDRESS:-127.0.0.1}\"",
            "UI_PORT=\"${ONT_UI_PORT:-8502}\"",
            'SERVER_ROOT="${ONT_SERVER_ROOT:-/data}"',
            'SERVER_START="${ONT_SERVER_START:-$SERVER_ROOT}"',
            'LOG_ROOT="${ONT_LOG_ROOT:-$PROJECT_DIR/usage_logs}"',
            'UI_RUN_ROOT="${ONT_UI_RUN_ROOT:-$PROJECT_DIR/ui_runs}"',
            "./launch_ui.sh --restart",
            '"$LOG_ROOT/ui_access.log"',
        ):
            self.assertIn(expected, text)

    def test_mcet13_launcher_uses_public_pc_paths(self):
        text = (ROOT / "ONT_Plasmid_Analyzer_MCET13_One_Click.bat").read_text(
            encoding="utf-8"
        )
        for expected in (
            'set "SERVER_USER=MCET13"',
            'set "REMOTE_PORT=8513"',
            'set "REMOTE_PROJECT=/home/MCET13/NGS_ONT"',
            'set "REMOTE_ENV=/home/MCET13/conda_envs/NGS_ONT_env"',
            'set "SERVER_ROOT=/data/user"',
            'set "SERVER_START=/data/user/MCET13"',
            'set "LOG_ROOT=/home/MCET13/.ont_plasmid_analyzer/logs"',
            'set "UI_RUN_ROOT=/data/user/MCET13/ui_runs"',
            "ONT_PROJECT_DIR=%REMOTE_PROJECT%",
            "ONT_CONDA_ENV=%REMOTE_ENV%",
            "ONT_SERVER_ROOT=%SERVER_ROOT%",
            "ONT_SERVER_START=%SERVER_START%",
            "ONT_LOG_ROOT=%LOG_ROOT%",
            "ONT_UI_RUN_ROOT=%UI_RUN_ROOT%",
        ):
            self.assertIn(expected, text)

    def test_pipeline_avoids_legacy_ont_profile_that_skips_indels(self):
        text = (ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        self.assertIn('BCFTOOLS_PLATFORM_MODE="ont-sup"', text)
        self.assertNotIn('BCFTOOLS_PLATFORM_OPT="-X ont"', text)
        self.assertIn('bcftools mpileup "${BCFTOOLS_PLATFORM_ARGS[@]}" -Ou', text)

    def test_launch_ui_forwards_server_browser_and_storage_paths(self):
        text = (ROOT / "launch_ui.sh").read_text(encoding="utf-8")
        for expected in (
            'SERVER_START="${ONT_SERVER_START:-$SERVER_ROOT}"',
            'ONT_SERVER_ROOT="$SERVER_ROOT" ONT_SERVER_START="$SERVER_START"',
            'ONT_LOG_ROOT="$LOG_ROOT" ONT_UI_RUN_ROOT="$UI_RUN_ROOT"',
        ):
            self.assertIn(expected, text)

    def test_pipeline_exports_repository_path_to_parallel_workers(self):
        text = (ROOT / "run_pipeline.sh").read_text(encoding="utf-8")
        self.assertIn("export SCRIPT_DIR DATA_ROOT", text)
        self.assertIn('$SCRIPT_DIR/ont_ui/fastq_qc.py', text)


if __name__ == "__main__":
    unittest.main()
