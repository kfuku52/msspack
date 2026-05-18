from __future__ import annotations

from .utils import MSSPackError


class ConfigError(MSSPackError):
    """Raised when a config file is invalid or missing required keys."""
