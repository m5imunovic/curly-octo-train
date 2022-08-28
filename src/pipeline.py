from itertools import chain

import hydra
from omegaconf import DictConfig


def run_graph_step(cfg: DictConfig):
    if cfg is not None:
        from graph.db_graph import run as run_db_graph
        run_db_graph(cfg)


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
    if 'reference' in cfg:
        run_reference_step(cfg.reference)

    if 'reads' in cfg:
        if 'reference' in cfg and cfg.reference is not None:
            if 'name' in cfg.reference:
                cfg.reads.species = cfg.reference.name
        run_generate_reads_step(cfg.reads)

    if 'asm' in cfg:
        if 'reference' in cfg and cfg.reference is not None:
            if 'name' in cfg.reference:
                cfg.asm.species = cfg.reference.name
            if 'chromosomes' in cfg.reference:
                import utils.path_helpers as ph
                from asm.assembler import assembly_experiment_path
                from reference.genome_generator import ref_chromosomes_path

                # TODO: check if hydra sweeper can be used to run assembly on multiple chromosomes
                chr_names = list(cfg.reference.chromosomes)
                chr_paths = [ref_chromosomes_path(cfg.reference) / f'{chr_name}.fasta' for chr_name in chr_names]
                for chr_path in chr_paths:
                    cfg.asm.params.long.reference = chr_path
                    reads_path = ph.get_simulated_data_path() / cfg.reference.name / chr_path.stem
                    out_path = assembly_experiment_path(cfg.asm) / chr_path.stem
                    run_assembly_step(cfg.asm, **{'reads_path': reads_path, 'out_path': out_path})
        else:
            run_assembly_step(cfg.asm)

    if 'graph' in cfg:
        if 'asm' in cfg and cfg.asm is not None:
            from asm.assembler import assembly_experiment_path
            asm_experiment_root = assembly_experiment_path(cfg.asm)
            chr_paths = list(filter(lambda x: x.is_dir(), asm_experiment_root.iterdir()))
            exp_paths = [list(filter(lambda x: x.is_dir(), chr_path.iterdir())) for chr_path in chr_paths]
            for exp_path in chain(*exp_paths):
                cfg.graph.gfa_path = exp_path / 'graph.gfa'
                cfg.graph.mult_info_path = exp_path / 'mult.info'
                run_graph_step(cfg.graph)
        else:
            run_graph_step(cfg.graph)


@hydra.main(version_base=None, config_path='../config', config_name='config')
def main(cfg: DictConfig):
    run(cfg)


if __name__ == '__main__':
    main()
