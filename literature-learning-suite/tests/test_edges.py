import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import gen_edges  # noqa: E402


class EdgeTests(unittest.TestCase):
    def setUp(self):
        self.papers = [
            {
                "id": "A",
                "entities": {"genes": ["TP53"], "methods": ["RNA-seq"]},
                "analysis": {
                    "cross_references": [
                        {
                            "ref_id": "B",
                            "relation": "extends",
                            "description": "Adds an independent model.",
                        }
                    ]
                },
            },
            {
                "id": "B",
                "entities": {"genes": ["TP53"], "methods": ["RNA-seq"]},
            },
        ]

    def test_cross_references_extraction(self):
        """Tier7 cross_references should be extracted from paper records."""
        refs = list(gen_edges.cross_references(self.papers[0]))
        self.assertGreaterEqual(len(refs), 1)
        self.assertEqual(refs[0]["ref_id"], "B")

    def test_curated_entities_extraction(self):
        """Curated entities should include genes and methods."""
        entities = gen_edges.curated_entities(self.papers[0])
        self.assertIn("TP53", str(entities))

    def test_analysis_block_extraction(self):
        """Analysis field should be accessible."""
        block = gen_edges.analysis_block(self.papers[0])
        self.assertIsInstance(block, dict)

    def test_entity_block_extraction(self):
        """Entity block should be accessible."""
        block = gen_edges.entity_block(self.papers[0])
        self.assertIsInstance(block, dict)


if __name__ == "__main__":
    unittest.main()
