import hydra
from omegaconf import DictConfig

import utils.path_helpers as ph
from reference.genome_generator import run


@hydra.main(version_base=None, config_path=str(ph.get_config_root()), config_name="reference")
def main(cfg: DictConfig):
    ph.adjust_cfg_paths(cfg)
    run(cfg)


if __name__ == "__main__":
    main()
