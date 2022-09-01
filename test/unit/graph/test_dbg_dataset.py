from cgi import test
from graph.dbg_dataset import DBGDataset

from torch_geometric.data import Data


def test_use_dbgdataset(test_datasets_root):
    dataset_root = test_datasets_root / 'random_species_09_01'
    dataset = DBGDataset(root=dataset_root)
    assert len(dataset) == 2
    ds0 = dataset[0]
    assert type(ds0) == Data
    assert ds0.num_nodes == 808
    assert ds0.num_edges == 644
    assert dataset is not None