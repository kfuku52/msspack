from __future__ import annotations

from pathlib import Path

from .config_errors import ConfigError
from .config_loading import load_sections, read_config_data
from .config_models import (
    BuscoConfig,
    InputsConfig,
    MSSPackConfig,
    PipelineConfig,
    ProjectConfig,
    ReferenceConfig,
    SampleConfig,
    StCommentConfig,
    SubmissionConfig,
    SubmitterConfig,
    ToolsConfig,
)
from .config_validation import validate_config

__all__ = [
    "ConfigError",
    "BuscoConfig",
    "InputsConfig",
    "MSSPackConfig",
    "PipelineConfig",
    "ProjectConfig",
    "ReferenceConfig",
    "SampleConfig",
    "StCommentConfig",
    "SubmissionConfig",
    "SubmitterConfig",
    "ToolsConfig",
    "load_config",
]


def load_config(path: str | Path) -> MSSPackConfig:
    config_path = Path(path).expanduser().resolve()
    config = MSSPackConfig(
        base_dir=config_path.parent,
        **load_sections(read_config_data(config_path)),
    )
    validate_config(config)
    return config
