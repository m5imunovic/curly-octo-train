from unittest import mock

from hydra import compose, initialize_config_dir
from omegaconf import OmegaConf

from pipeline import run as run_pipeline
from utils.path_helpers import get_default_cfg_path


@mock.patch('pipeline.run_graph_step')
@mock.patch('pipeline.run_assembly_step')
@mock.patch('pipeline.run_generate_reads_step')
@mock.patch('pipeline.run_reference_step')
def test_run_pipeline(run_reference_step, run_generate_reads_step, run_assembly_step, run_graph_step):
    cfg  = OmegaConf.create({
        'reference': None,
        'reads': None,
        'asm': None,
        'graph': None
    })


    run_pipeline(cfg)

    run_reference_step.assert_called_once()
    run_generate_reads_step.assert_called_once()
    run_assembly_step.assert_called_once()
    run_graph_step.assert_called_once()


@mock.patch('pipeline.run_graph_step')
@mock.patch('pipeline.run_assembly_step')
@mock.patch('pipeline.run_generate_reads_step')
@mock.patch('pipeline.run_reference_step')
def test_skipped_reference_and_reads_steps(run_reference_step, run_generate_reads_step, run_assembly_step, run_graph_step):
    cfg  = OmegaConf.create({
        'asm': None
    })


    run_pipeline(cfg)

    run_reference_step.assert_not_called()
    run_generate_reads_step.assert_not_called()
    run_assembly_step.assert_called_once()
    run_graph_step.assert_not_called()


@mock.patch('pipeline.run_graph_step')
@mock.patch('pipeline.run_assembly_step')
@mock.patch('pipeline.run_generate_reads_step')
@mock.patch('pipeline.run_reference_step')
def test_replaced_config_values(run_reference_step, run_generate_reads_step, run_assembly_step, run_graph_step):
    # TODO: replace default config with test pipeline config
    cfg_path = get_default_cfg_path()
    with initialize_config_dir(str(cfg_path.parent), version_base=None):
        cfg = compose(str(cfg_path.name), overrides=['reference.name=dummy', 'asm.experiment=dummy_experiment'])
        run_pipeline(cfg)
        reference_cfg = run_reference_step.call_args.args[0]
        assert reference_cfg.name == 'dummy'
        simulator_cfg = run_generate_reads_step.call_args.args[0]
        assert simulator_cfg.species == 'dummy'
        asm_cfg = run_assembly_step.call_args.args[0]
        assert asm_cfg.species == 'dummy'
        assert 'dummy/chromosomes/' in str(asm_cfg.params.long.reference)
        graph_cfg = run_graph_step.call_args.args[0]
        assert 'dummy_experiment/graph.gfa' in str(graph_cfg.gfa_path)
