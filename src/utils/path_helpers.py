from pathlib import Path


def get_project_root() -> Path:
    return Path(__file__).parent.parent.parent


def get_vendor_path() -> Path:
    return get_project_root() / 'vendor'


def get_data_path() -> Path:
    return get_project_root() / 'data'


def get_ref_path() -> Path:
    return get_data_path() / 'references'


def get_simulated_data_path() -> Path:
    return get_data_path() / 'simulated'


def get_assemblies_path() -> Path:
    return get_data_path() / 'assemblies'


def get_default_cfg_path() -> Path:
    return get_project_root() / 'config' / 'config.yaml'


def project_root_append(path: str):
    return get_project_root() / path

