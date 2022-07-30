from unittest import mock
from omegaconf import OmegaConf

from pipeline import run as run_pipeline


@mock.patch('pipeline.run_assembly_step')
@mock.patch('pipeline.run_generate_reads_step')
@mock.patch('pipeline.run_reference_step')
def test_run_pipeline(run_reference_step, run_generate_reads_step, run_assembly_step):
    cfg  = OmegaConf.create({
        'reference': None,
        'reads': None,
        'asm': None
    })


    run_pipeline(cfg)

    run_reference_step.assert_called_once()
    run_generate_reads_step.assert_called_once()
    run_assembly_step.assert_called_once()


@mock.patch('pipeline.run_assembly_step')
@mock.patch('pipeline.run_generate_reads_step')
@mock.patch('pipeline.run_reference_step')
def test_skipped_reference_and_reads_steps(run_reference_step, run_generate_reads_step, run_assembly_step):
    cfg  = OmegaConf.create({
        'asm': None
    })


    run_pipeline(cfg)

    run_reference_step.assert_not_called()
    run_generate_reads_step.assert_not_called()
    run_assembly_step.assert_called_once()