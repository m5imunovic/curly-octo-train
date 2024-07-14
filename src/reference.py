import hydra
from omegaconf import DictConfig

import utils.path_helpers as ph
from reference.genome_generator import run
from reference.bed_genome_generator import run as generate_bed_reference


@hydra.main(version_base=None, config_path=str(ph.get_config_root()), config_name="reference")
def main(cfg: DictConfig):
    ph.adjust_cfg_paths(cfg)
    run(cfg)
    if "bed" in cfg.reference and cfg.reference.bed is not None:
        generate_bed_reference(cfg)


if __name__ == "__main__":
    main()
