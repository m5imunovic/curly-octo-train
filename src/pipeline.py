
import hydra
from omegaconf import DictConfig


def run_graph_step(cfg: DictConfig):
    if cfg is not None:
        from graph.db_graph import run as run_db_graph
        run_db_graph(cfg)


def run_assembly_step(cfg: DictConfig):
    if cfg is not None:
        from asm.assembler import run as run_assembly
        run_assembly(cfg)


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
                from reference.genome_generator import ref_chromosomes_path
                # TODO: introduce sweeper to run assembly on multiple chromosomes
                chr_name = list(cfg.reference.chromosomes)[0]
                chr_path = ref_chromosomes_path(cfg.reference) / f'{chr_name}.fasta'
                cfg.asm.params.long.reference = chr_path

        run_assembly_step(cfg.asm)

    if 'graph' in cfg:
        if 'asm' in cfg and cfg.asm is not None:
            from asm.assembler import assembly_experiment_path
            cfg.graph.gfa_path = assembly_experiment_path(cfg.asm) / 'graph.gfa'
        run_graph_step(cfg.graph)


@hydra.main(version_base=None, config_path='../config', config_name='config')
def main(cfg: DictConfig):
    run(cfg)


if __name__ == '__main__':
    main()