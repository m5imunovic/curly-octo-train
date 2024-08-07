from pathlib import Path

from omegaconf import DictConfig, OmegaConf

# TODO: clean up path helpers as we don't need all of them, also, make config root configurable with envvar


def get_project_root() -> Path:
    return Path(__file__).parent.parent.parent


def get_vendor_path() -> Path:
    return get_project_root() / "vendor"


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
    for key, value in cfg.paths.items():
        if isinstance(value, str):
            cfg.paths[key] = project_root / value
        elif isinstance(value, list):
            cfg.paths[key] = [project_root / path for path in value]
        else:
            raise ValueError(f"Unknown type for path {key}: {type(value)}")
