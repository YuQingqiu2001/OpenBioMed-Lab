import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import gen_edges  # noqa: E402
import selfcheck_knowledge_graph as selfcheck  # noqa: E402
import build_network as network  # noqa: E402


class HermesCompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.nested = {
            "id": "PMID:1",
            "title": "Nested record",
            "analysis": {
                "core_question": "How does TP53 alter the response?",
                "subquestions": ["q1", "q2", "q3", "q4", "q5"],
                "ces_chains": [
                    {
                        "claim": f"claim {index}",
                        "evidence": "specific evidence with enough detail for audit",
                        "synthesis": "synthesis",
                    }
                    for index in range(5)
                ],
                "mechanism": {"steps": ["TP53", "DNA damage", "apoptosis"]},
                "hidden_axes": [{"observation": "obs", "interpretation": "axis"}],
                "conceptual_contribution": {"new_concepts": ["concept"]},
                "cross_references": [
                    {
                        "ref_id": "PMID:2",
                        "relation": "extends",
                        "description": "Tests the same mechanism in another model.",
                    }
                    for _ in range(5)
                ],
            },
            "entities": {
                "genes": ["TP53"],
                "pathways": ["apoptosis"],
                "diseases": ["cancer"],
                "methods": ["single-cell"],
            },
        }
        self.legacy = {
            "id": "PMID:2",
            "title": "Legacy record",
            "tier2_core_question": "What is the legacy question?",
            "tier2_subquestions": ["q1"],
            "tier3_ces_chains": [{"evidence": "legacy evidence"}],
            "tier4_mechanism_cascade": {"cascade": ["A", "B"]},
            "tier5_hidden_axis": [{"observation": "legacy"}],
            "tier6_concept_innovation": {"new_concepts": ["legacy concept"]},
            "tier7_cross_refs": [{"ref_id": "PMID:1", "relation": "supports"}],
            "genes": ["EGFR"],
            "technologies": ["spatial transcriptomics"],
        }

    def test_nested_record_maps_to_hermes_runtime(self):
        self.assertEqual(
            gen_edges.core_question(self.nested),
            "How does TP53 alter the response?",
        )
        self.assertEqual(
            gen_edges.mechanism_steps(self.nested),
            ["TP53", "DNA damage", "apoptosis"],
        )
        entities = gen_edges.curated_entities(self.nested)
        self.assertIn("TP53", entities["molecules"])
        self.assertIn("single-cell", entities["methods"])
        self.assertIn("specific evidence", gen_edges.paper_text(self.nested))

    def test_legacy_record_remains_supported(self):
        self.assertEqual(
            gen_edges.core_question(self.legacy),
            "What is the legacy question?",
        )
        self.assertEqual(gen_edges.mechanism_steps(self.legacy), ["A", "B"])
        entities = gen_edges.curated_entities(self.legacy)
        self.assertIn("EGFR", entities["molecules"])
        self.assertIn("spatial transcriptomics", entities["methods"])

    def test_selfcheck_audits_nested_s_tier(self):
        self.assertEqual(selfcheck.analysis_tier(self.nested), "S")
        view = selfcheck.analysis_view(self.nested)
        self.assertEqual(len(view["ces_chains"]), 5)
        self.assertTrue(selfcheck.has_content(view["mechanism"]))

    def test_network_reads_nested_and_legacy_findings(self):
        nested_claims = network.paper_claims(self.nested)
        legacy_claims = network.paper_claims(self.legacy)
        self.assertEqual(len(nested_claims), 5)
        self.assertEqual(len(legacy_claims), 1)
        self.assertIn(
            "How does TP53",
            network.paper_finding(self.nested, nested_claims),
        )
        self.assertIn(
            "What is the legacy",
            network.paper_finding(self.legacy, legacy_claims),
        )


if __name__ == "__main__":
    unittest.main()
