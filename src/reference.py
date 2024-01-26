import hydra
from omegaconf import OmegaConf

import utils.path_helpers as ph
from utils.config_cs import Config, register_configs

register_configs()


@hydra.main(version_base=None, config_path=str(ph.get_config_root()), config_name="reference")
def main(cfg: Config):
    ph.adjust_cfg_paths(cfg)
    print(OmegaConf.to_yaml(cfg))
    # run(cfg)


if __name__ == "__main__":
    main()
