from pathlib import Path

import pytorch_lightning as pl
import torch
import torch.nn as nn
import torch_geometric.transforms as T
from torch_geometric.loader import DataLoader

from graph.dbg_dataset import DBGDataset
from models.model_di_graph import ModuleDiGraph


LR = 1e-3
#DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
DEVICE = torch.device('cpu')
NODE_FEATURES=6
HIDDEN_NODE_FEATURES=8
EPOCHS = 200
PATIENCE = 2
DECAY = 0.95


class DBGLightningModule(pl.LightningModule):
    def __init__(self, model: nn.Module, pos_weight: torch.Tensor):
        super().__init__()
        self.model = model
        self.criterion = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        self.batch_size = 1

    def configure_optimizers(self) -> torch.optim.Optimizer:
        optimizer = torch.optim.Adam(self.parameters(), lr=LR)
        lr_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=DECAY, patience=PATIENCE, verbose=True)
        return { "optimizer": optimizer, "lr_scheduler": lr_scheduler, "monitor": "val_loss" }

    def training_step(self, batch, batch_idx):
        scores = self.model(batch.x.float(), batch.edge_index)
        expected_scores = batch.y.float().unsqueeze(-1)
        loss = self.criterion(scores, expected_scores)

        self.log("train_loss", loss, on_step=True, on_epoch=True, prog_bar=True, logger=True, batch_size=self.batch_size)
        return loss

    def validation_step(self, batch, batch_idx):
        scores = self.model(batch.x.float(), batch.edge_index)
        expected_scores = batch.y.float().unsqueeze(-1)
        loss = self.criterion(scores, expected_scores)

        self.log("val_loss", loss, on_step=True, on_epoch=True, prog_bar=True, logger=True, batch_size=self.batch_size)
        return loss

    def test_step(self, batch, batch_idx):
        scores = model(batch.x.float(), batch.edge_index)
        expected_scores = batch.y.float().unsqueeze(-1)
        scores_rounded = torch.round(torch.sigmoid(scores))
        correct = (scores_rounded == expected_scores).sum().item()
        test_acc = int(correct) / int(len(batch.y))

        self.log("test_acc", test_acc, on_step=True, on_epoch=True, prog_bar=True, logger=True, batch_size=self.batch_size)
        return test_acc

if __name__ == '__main__':

    path_train = (Path(__file__).parent / '../../data/datasets/random_species_09_14').resolve()
    path_validate = (Path(__file__).parent / '../../data/datasets/random_species_09_15').resolve()
    path_test = (Path(__file__).parent / '../../data/datasets/random_species_09_17').resolve()
    train_ds = DBGDataset(root=path_train, transform=T.NormalizeFeatures(attrs=['x']))
    validate_ds = DBGDataset(root=path_validate, transform=T.NormalizeFeatures(attrs=['x']))
    test_ds = DBGDataset(root=path_test, transform=T.NormalizeFeatures(attrs=['x']))
    pos_to_neg_ratio = torch.mean(torch.tensor([(g.y == 1).sum() / (g.y == 0).sum() for g in train_ds])).item()
    pos_weight = torch.tensor([1/pos_to_neg_ratio], device=DEVICE)

    model = ModuleDiGraph(node_features=NODE_FEATURES, hidden_features=HIDDEN_NODE_FEATURES, device=DEVICE)

    pl_model = DBGLightningModule(model=model, pos_weight=pos_weight)

    trainer = pl.Trainer(
        check_val_every_n_epoch=1,
        max_epochs=EPOCHS, accelerator='gpu',
        devices=1,
        log_every_n_steps=1,
    )
    trainer.fit(
        model=pl_model,
        train_dataloaders=DataLoader(train_ds, num_workers=4, batch_size=1),
        val_dataloaders=DataLoader(validate_ds, num_workers=4, batch_size=1)
    )

    trainer.test(
        model=pl_model,
        dataloaders=DataLoader(test_ds, num_workers=4, batch_size=1)
    )