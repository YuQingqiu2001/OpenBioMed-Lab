import tempfile
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from journal_metrics import import_metrics, lookup, normalize_record  # noqa: E402


class JournalMetricTests(unittest.TestCase):
    def test_import_and_lookup_2024(self):
        source = ROOT / "tests" / "fixtures" / "journal_metrics_2024.synthetic.json"
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder)
            metadata = import_metrics(source, target, 2024, "synthetic fixture")
            self.assertEqual(metadata["metric_year"], 2024)
            self.assertEqual(metadata["record_count"], 1)
            matches = lookup(target, "0000-0001", 5)
            self.assertEqual(matches[0]["metric_year"], 2024)
            self.assertEqual(matches[0]["jif"], 4.2)

    def test_case_insensitive_field_names(self):
        record = normalize_record({"Journal": "Example", "JIF": "3.1"}, 2024)
        self.assertEqual(record["name"], "Example")
        self.assertEqual(record["jif"], 3.1)


if __name__ == "__main__":
    unittest.main()
