import hydra
from omegaconf import DictConfig

from graph.construct_graph import construct_graph


def run(cfg: DictConfig):
    g = construct_graph(cfg)
    print(f'Number of edges {g.number_of_edges()}')
    print(f'Number of nodes {g.number_of_nodes()}')


@hydra.main(version_base=None, config_path='../../config/graph', config_name='db_graph')
def main(cfg: DictConfig):
    run(cfg)


if __name__ == '__main__':
    main()