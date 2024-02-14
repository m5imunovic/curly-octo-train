import hydra
from omegaconf import DictConfig, OmegaConf

import utils.path_helpers as ph
from reads.sample_reads_fastq import run


@hydra.main(version_base=None, config_path=str(ph.get_config_root()), config_name="sprofile")
def main(cfg: DictConfig):
    ph.adjust_cfg_paths(cfg)
    OmegaConf.resolve(cfg)
    run(cfg)


if __name__ == "__main__":
    main()
