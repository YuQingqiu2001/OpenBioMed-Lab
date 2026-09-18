from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_two_model_families_are_present():
    assert (ROOT / "models/cell_topology/weights/cell_topology_model.npz").is_file()
    assert all((ROOT / f"models/patch_feature/weights/fold_{fold}/survival_aggregator.pt").is_file() for fold in range(5))


def test_third_party_assets_are_not_distributed():
    assert not (ROOT / "weights").exists()
    assert not (ROOT / "vendor").exists()
    assert not list(ROOT.rglob("*.pth"))


def test_official_upstream_links_are_documented():
    notices = (ROOT / "docs/THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    assert "https://github.com/TIO-IKIM/CellViT-plus-plus" in notices
    assert "https://github.com/mahmoodlab/UNI" in notices
