from typing import Optional

import torch
import torch.nn as nn
from torch_geometric.nn.conv import ResGatedGraphConv


class ModuleDiGraph(nn.Module):
    def __init__(self,
                 node_features: int,
                 hidden_features: int,
                 device: Optional[torch.device] = None,
                 batch_norm: bool = False):
        super().__init__()

        self.W1 = nn.Linear(node_features, hidden_features, bias=True, device=device)
        self.W2 = nn.Linear(hidden_features, hidden_features, bias=True, device=device)

        self.gate = ResGatedGraphConv(hidden_features, hidden_features)

        self.scorer = nn.Linear(hidden_features, out_features=1, bias=True, device=device)

        self.reset_parameters()


    def reset_parameters(self):
        self.W1.reset_parameters()
        self.W2.reset_parameters()
        self.scorer.reset_parameters()

    def forward(self, x, edge_index):
        h = self.W1(x)
        h = torch.relu(h)
        h = self.W2(h)
        h = self.gate(h, edge_index)
        score = self.scorer(h)

        return score


if __name__ == '__main__':
    from pathlib import Path
    from graph.dbg_dataset import DBGDataset
    # from torch_geometric.transforms import NormalizeFeatures
    import torch_geometric.transforms as T

    path = (Path(__file__).parent / '../../data/datasets/random_species_09_10').resolve()
    dataset = DBGDataset(root=path, transform=T.NormalizeFeatures(attrs=['x']))

    sigmoid = nn.Sigmoid()

    #pos_to_neg_ratio = sum([((g.edata['y']==1).sum() / (g.edata['y']==0).sum()).item() for idx, g in ds_train]) / len(ds_train)

    pos_to_neg_ratio = torch.mean(torch.tensor([(g.y == 1).sum() / (g.y == 0).sum() for g in dataset])).item()

    for entry in dataset:
        data = entry
        #module = ResGatedGraphConv(6, 1)
        module = ModuleDiGraph(node_features=6, hidden_features=64)
        module.reset_parameters()
        scores = module.forward(data.x.float(), data.edge_index)
        y_hat = sigmoid(scores)
        print(y_hat.shape)
