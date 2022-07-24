import utils.path_helpers as ph


def test_get_project_root():
    root_path = ph.get_project_root()
    assert root_path.exists()
    toplevel_readme_file = root_path / "README.md"
    assert toplevel_readme_file.exists()

def test_get_vendor_path():
    vendor_path = ph.get_vendor_path()
    assert vendor_path.exists()
    assert vendor_path.name == "vendor"


def test_get_references_path():
    ref_path = ph.get_ref_path()
    assert ref_path.exists()
    assert ref_path.name == "references"


def test_get_default_cfg_path():
    cfg_path = ph.get_default_cfg_path()
    assert cfg_path.exists()
    assert cfg_path.name == "config.yaml"
