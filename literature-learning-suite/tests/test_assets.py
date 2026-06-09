import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from check_assets import check_assets  # noqa: E402
from init_workspace import seed_workspace_data  # noqa: E402
from validate_records import ALLOWED_RELATIONS, validate_edge  # noqa: E402


class AssetTests(unittest.TestCase):
    def test_bundled_assets_match_manifest(self):
        report = check_assets()
        self.assertEqual(report["failures"], [])
        counts = {item["path"]: item["records"] for item in report["datasets"]}
        self.assertEqual(counts["assets/data/bioc_genes.json"], 90125)
        self.assertEqual(counts["assets/data/kegg_pathways.json"], 25939)

    def test_runtime_relations_are_declared(self):
        expected = {
            "cites",
            "defines_concept",
            "shares_curated_entities",
            "shares_disease_method",
            "shares_molecules",
            "shares_paradigm",
            "shares_topic",
        }
        self.assertTrue(expected.issubset(ALLOWED_RELATIONS))
        self.assertIn(
            "forbidden relation",
            validate_edge(
                {
                    "source": "A",
                    "target": "B",
                    "relation": "same_journal",
                    "description": "Not a scientific relation.",
                }
            ),
        )

    def test_workspace_seed_preserves_local_data_until_refresh(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            copied = seed_workspace_data(root)
            self.assertIn("bioc_genes.json", copied)
            genes = root / "data" / "bioc_genes.json"
            self.assertEqual(len(json.loads(genes.read_text(encoding="utf-8"))), 90125)

            genes.write_text('["LOCAL"]\n', encoding="utf-8")
            self.assertNotIn("bioc_genes.json", seed_workspace_data(root))
            self.assertEqual(json.loads(genes.read_text(encoding="utf-8")), ["LOCAL"])

            self.assertIn("bioc_genes.json", seed_workspace_data(root, refresh=True))
            self.assertEqual(len(json.loads(genes.read_text(encoding="utf-8"))), 90125)


if __name__ == "__main__":
    unittest.main()
