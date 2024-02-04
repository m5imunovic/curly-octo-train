import hydra
from omegaconf import DictConfig

import utils.path_helpers as ph
from experiment.schedule import run


@hydra.main(version_base=None, config_path=str(ph.get_config_root()), config_name="experiment")
def main(cfg: DictConfig):
    ph.adjust_cfg_paths(cfg)
    run(cfg)


if __name__ == "__main__":
    main()
