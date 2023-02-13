import hydra
from omegaconf import DictConfig

import utils.path_helpers as ph


def run_graph_step(cfg: DictConfig, **kwargs):
    if cfg is not None:
        from graph.db_graph import run as run_db_graph

        run_db_graph(cfg, **kwargs)


def run_assembly_step(cfg: DictConfig, **kwargs):
    if cfg is not None:
        from asm.assembler import run as run_assembly

        run_assembly(cfg, **kwargs)


def run_generate_reads_step(cfg: DictConfig):
    if cfg is not None:
        from reads.reads_simulator import run as run_generate_reads

        run_generate_reads(cfg)


def run_reference_step(cfg: DictConfig):
    if cfg is not None:
        from reference.genome_generator import run as run_reference_step

        run_reference_step(cfg)


def run(cfg: DictConfig):
    if "reference" in cfg:
        run_reference_step(cfg)

    if "reads" in cfg:
        run_generate_reads_step(cfg)

    if "asm" in cfg:
        if "reference" in cfg and cfg.reference is not None:
            # TODO: overwrite should be inspected in run command to see if anything should be skipped
            if "chromosomes" in cfg.reference and cfg.asm.overwrite:
                from asm.assembler import assembly_experiment_path
                from reference.genome_generator import ref_chromosomes_path

                # TODO: check if hydra sweeper can be used to run assembly on multiple chromosomes
                # TODO: Code should infer the chromosomes from the reference genome path and not from the config
                chr_names = list(cfg.reference.chromosomes)
                chr_paths = [ref_chromosomes_path(cfg) / f"{chr_name}.fasta" for chr_name in chr_names]
                for chr_path in chr_paths:
                    cfg.asm.params.long.reference = chr_path
                    reads_path = cfg.paths.reads_dir / cfg.species_name / chr_path.stem
                    out_path = assembly_experiment_path(cfg) / chr_path.stem
                    run_assembly_step(cfg, **{"reads_path": reads_path, "out_path": out_path})
        else:
            run_assembly_step(cfg)

    if "graph" in cfg:
        if "asm" in cfg and cfg.asm is not None:
            cfg.graph.experiment = cfg.asm.experiment
            from asm.assembler import assembly_experiment_path

            asm_experiment_root = assembly_experiment_path(cfg)
            out_path = cfg.paths.datasets_dir / cfg.graph.experiment
            run_graph_step(cfg, **{"assemblies_path": asm_experiment_root, "out_path": out_path})
        else:
            run_graph_step(cfg)


@hydra.main(version_base="1.2", config_path="../config", config_name="config")
def main(cfg: DictConfig):
    ph.adjust_cfg_paths(cfg)
    run(cfg)


if __name__ == "__main__":
    main()
