import hydra
from omegaconf import DictConfig

import utils.path_helpers as ph
from graph.db_graph import run


@hydra.main(version_base=None, config_path=str(ph.get_config_root()), config_name="graph.yaml")
def main(cfg: DictConfig):
    ph.adjust_cfg_paths(cfg)
    run(cfg)


if __name__ == "__main__":
    main()
