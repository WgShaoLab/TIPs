import yaml
from dataclasses import dataclass
from typing import Any, Dict



@dataclass
class PipelineConfig:
    """
    Lightweight wrapper for YAML config.
    You can add validation logic here later if needed.
    """
    raw: Dict[str, Any]

    @property
    def sample_path(self) -> str:
        return self.raw["sample"]["path"]

    @property
    def gpu_id(self) -> str:
        return self.raw["gpu"]["id"]

    # You can keep adding convenience properties for commonly-used fields.
    # For less common values, just access config.raw[...] directly.


def load_config(path: str) -> PipelineConfig:
    """
    Load user YAML config into a PipelineConfig instance.
    """
    with open(path, "r") as f:
        data = yaml.safe_load(f)
    return PipelineConfig(raw=data)
