import hydra
import networkx as nx
from omegaconf import OmegaConf
from typeguard import typechecked

from graph.load_gfa import load_gfa, GFAFeatures
from graph.load_mult_info import load_mult_info, compare_edge_ids
from graph.gfa_parser import parse_gfa


class DBGraph(nx.DiGraph):
    @staticmethod
    @typechecked
    def from_gfa_features(gfa_features: GFAFeatures) -> nx.DiGraph:
        g = nx.DiGraph()
        for src, tgt_dict in gfa_features['edges'].items():
            for tgt, features in tgt_dict.items():
                g.add_edge(src, tgt, **features)

        g.add_nodes_from(gfa_features['nodes'].items())

        return nx.line_graph(g)


def run(cfg: OmegaConf):
    ref_path = cfg.ref_path
    gfa_path = cfg.gfa_path

    gfa_features = parse_gfa(gfa_path)
    gfa_features = load_gfa(gfa_path, skip_rc=True)
    mult_info_path = cfg.mult_info_path
    mult_info = load_mult_info(mult_info_path)
    ##g = DBGraph.from_gfa_features(gfa_features)
    is_diff_left, is_diff_right = compare_edge_ids(mult_info.keys(), gfa_features['nodes'].keys())
    print(f'{is_diff_left=} {is_diff_right=}')
    #print(gfa_features)

#    New chromosome chr1(200000)                                                                                                                                                                            
#    [0, 4624] -> 2586835047048751220864862834605086160691T [0, 4624]
#    [4624, 9065] -> 1566937966847838364864340278423830827451T [0, 4441]
#    [9065, 9267] -> 16535373539295126568437052775843160720G [0, 202]
#    [9267, 13105] -> 2777391822930213979519031551135732399490T [0, 3838]
#    [13105, 13199] -> 736441387725877665520783138001251207240A [0, 94]
#    [15481, 15505] -> 168642341276858936513619140475527997251A [0, 24]
#    [16534, 16769] -> 1600684769314781166715936728669035374080G [0, 235]


    
    from rolling_hash import RollingHash
    with open(ref_path, 'r') as f:
        for idx,line in enumerate(f):
            if idx == 0:
                continue
            line = line.strip()
            parts = [(0, 4624),
                     (4624, 9065), 
                     (9065, 9267), 
                     (9267, 13105), 
                     (13105, 13199),
                     (15481, 15505),
                     (16534, 16769)
            ]
                    
            rh = RollingHash()
            for part_start, part_end in parts:
                #print(line[part_start:part_end])
                print(rh.hash(line[part_start:part_end + max(part_end - part_start, 501)], 0))
                print(rh.hash(line[part_start:part_end + 501], 0))
            break

    # range: [0:4584]


@hydra.main(version_base=None, config_path='../../config/graph', config_name='db_graph')
def main(cfg: OmegaConf):
    run(cfg)


if __name__ == "__main__":
    main()