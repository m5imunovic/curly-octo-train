from utils import path_helpers as ph
from omegaconf import OmegaConf


OmegaConf.register_new_resolver('project_root', ph.project_root_append, replace=True)
