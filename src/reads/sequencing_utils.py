import logging
from pathlib import Path

import wandb
from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig, OmegaConf

logger = logging.getLogger(__name__)


def upload_pbsim3_profile_to_wandb(cfg: DictConfig) -> bool:
    """Uploads a pbsim3 profile to wandb.

    Args:
        profile: A pbsim3 profile.
        run_id: A wandb run id.
    """

    if cfg.wandb is None:
        logger.warn("WandB config is not defined, skipping artifact logging!")
        return False

    metadata = OmegaConf.to_container(cfg.wandb)
    profile_id_name = f"{cfg.reads.params.long['sample-profile-id']}.fastq"
    profile = Path(cfg.reads.profile.path) / f"sample_profile_{profile_id_name}"
    hydra_cfg_dir = Path(HydraConfig().get().run.dir).absolute() / ".hydra"

    logger.info(f"Uploading profile {profile} to wandb.")
    run = wandb.init(project="pbsim3", job_type="add-profile")
    artifact = wandb.Artifact(name=profile.stem, type="pbsim3-profile", incremental=True, metadata=metadata)
    artifact.add_file(local_path=profile)
    artifact.add_file(local_path=profile.with_suffix(".stats"))
    artifact.add_dir(local_path=str(hydra_cfg_dir), name="metadata")
    run.log_artifact(artifact)
    logger.info("Profile uploaded to wandb.")

    return True
