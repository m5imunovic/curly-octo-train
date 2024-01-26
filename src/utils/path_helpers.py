from pathlib import Path

from omegaconf import DictConfig, OmegaConf


def get_project_root() -> Path:
    return Path(__file__).parent.parent.parent


def get_vendor_path() -> Path:
    return get_project_root() / "vendor"


def get_data_path() -> Path:
    return get_project_root() / "data"


def get_ref_path() -> Path:
    return get_data_path() / "references"


def get_assemblies_path() -> Path:
    return get_data_path() / "assemblies"


def get_datasets_path() -> Path:
    return get_data_path() / "datasets"


def get_config_root() -> Path:
    return get_project_root() / "config"


def get_default_cfg_path() -> Path:
    return get_config_root() / "config.yaml"


def project_root_append(path: str):
    return get_project_root() / path


def adjust_cfg_paths(cfg: DictConfig):
    assert "paths" in cfg, 'Config must contain a "paths" section'

    project_root = get_project_root()
    OmegaConf.resolve(cfg.paths)
    for path_key, path_value in cfg.paths.items():
        if isinstance(path_value, Path):
            if not path_value.is_absolute():
                cfg.paths[path_key] = project_root / path_value
        else:
            raise ValueError(f"Unexpected type for path {path_key}: {type(path_value)}")
