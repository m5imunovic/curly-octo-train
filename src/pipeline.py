import hydra
from omegaconf import OmegaConf


def run_assembly_step(cfg: OmegaConf):
    if cfg is not None:
        from asm.assembler import run as run_assembly
        run_assembly(cfg)


def run_generate_reads_step(cfg: OmegaConf):
    if cfg is not None:
        from reads.reads_simulator import run as run_generate_reads
        run_generate_reads(cfg)


def run_reference_step(cfg: OmegaConf):
    if cfg is not None:
        from reference.genome_generator import run as run_reference_step
        run_reference_step(cfg)


def run(cfg: OmegaConf):
    if 'reference' in cfg:
        run_reference_step(cfg['reference'])

    if 'reads' in cfg:
        run_generate_reads_step(cfg['reads'])

    if 'asm' in cfg:
        run_assembly_step(cfg['asm'])



@hydra.main(version_base=None, config_path='../config', config_name='config')
def main(cfg: OmegaConf):
    run(cfg)


if __name__ == '__main__':
    main()