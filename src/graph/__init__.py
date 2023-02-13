from omegaconf import OmegaConf

from utils import path_helpers as ph

OmegaConf.register_new_resolver("project_root", ph.project_root_append, replace=True)
