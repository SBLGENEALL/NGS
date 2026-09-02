import csv
import tempfile
import unittest
from pathlib import Path

from ont_ui.audit import append_usage_event


class AuditLogTests(unittest.TestCase):
    def test_usage_events_append_with_project_user_and_timestamp(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "usage.csv"
            first = append_usage_event(
                path,
                event="STARTED",
                project_name="Project A",
                user_name="Jongin Baek",
                reference_count=5,
                analysis_count=15,
                job_id="job-1",
            )
            append_usage_event(
                path,
                event="COMPLETED",
                project_name="Project A",
                user_name="Jongin Baek",
                reference_count=5,
                analysis_count=15,
                job_id="job-1",
            )

            with path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["timestamp_kst"], first)
            self.assertEqual(rows[0]["project_name"], "Project A")
            self.assertEqual(rows[0]["user_name"], "Jongin Baek")
            self.assertEqual([row["event"] for row in rows], ["STARTED", "COMPLETED"])

    def test_spreadsheet_formula_prefix_is_neutralized(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "usage.csv"
            append_usage_event(
                path,
                event="FAILED",
                project_name="=unsafe",
                user_name="+unsafe",
                reference_count=0,
                analysis_count=0,
            )
            with path.open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["project_name"], "'=unsafe")
            self.assertEqual(row["user_name"], "'+unsafe")


if __name__ == "__main__":
    unittest.main()
