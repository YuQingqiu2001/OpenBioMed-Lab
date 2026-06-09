import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import literature_search  # noqa: E402


class SearchTests(unittest.TestCase):
    @patch("literature_search.request")
    def test_crossref_normalization(self, mocked_request):
        mocked_request.return_value = json.dumps(
            {
                "message": {
                    "items": [
                        {
                            "DOI": "10.1000/ABC",
                            "title": ["A test paper"],
                            "author": [{"given": "Ada", "family": "Lovelace"}],
                            "published": {"date-parts": [[2024, 1, 1]]},
                            "container-title": ["Test Journal"],
                            "URL": "https://doi.org/10.1000/ABC",
                            "type": "journal-article",
                        }
                    ]
                }
            }
        ).encode()
        records = literature_search.search_crossref("test", 1)
        self.assertEqual(records[0]["id"], "DOI:10.1000/abc")
        self.assertEqual(records[0]["authors"], ["Ada Lovelace"])
        self.assertEqual(records[0]["year"], 2024)

    @patch("literature_search.request")
    def test_preprint_date_window_and_filter(self, mocked_request):
        mocked_request.return_value = json.dumps(
            {
                "messages": [{"total": 1}],
                "collection": [
                    {
                        "doi": "10.1101/2024.01.01.000001",
                        "title": "Spatial transcriptomics example",
                        "abstract": "A reproducible test.",
                        "authors": "A Author; B Author",
                        "date": "2024-01-02",
                        "version": "1",
                        "category": "bioinformatics",
                    }
                ],
            }
        ).encode()
        records = literature_search.search_preprint(
            "biorxiv", "spatial transcriptomics", 5, "2024-01-01", "2024-01-31"
        )
        self.assertEqual(len(records), 1)
        self.assertIn("/2024-01-01/2024-01-31/0/json", mocked_request.call_args.args[0])


if __name__ == "__main__":
    unittest.main()
