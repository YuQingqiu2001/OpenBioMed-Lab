import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from ll_common import normalize_doi, normalize_title, stable_id  # noqa: E402


class CommonTests(unittest.TestCase):
    def test_normalize_doi(self):
        self.assertEqual(
            normalize_doi("https://doi.org/10.1000/Example."), "10.1000/example"
        )

    def test_stable_id_priority(self):
        self.assertEqual(stable_id({"pmid": "123"}), "PMID:123")
        self.assertEqual(stable_id({"doi": "10.1/ABC"}), "DOI:10.1/abc")

    def test_normalize_title(self):
        self.assertEqual(normalize_title("<b>A  Study!</b>"), "a study")


if __name__ == "__main__":
    unittest.main()
