The installation of Pytorch Geometric package is highly dependent on the system configuration.
Follow the instructions on official [website](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html) for
detailed instructions. The pip installation command will look something like this:

```bash
pip install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu113
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric -f https://data.pyg.org/whl/torch-1.12.1+cu113.html
```
