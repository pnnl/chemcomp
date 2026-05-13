"""
Configuration module for output settings and global parameters.
"""

import os
from dataclasses import dataclass
from typing import Optional
import warnings

# Suppress warnings globally
warnings.filterwarnings('ignore')


@dataclass
class OutputConfig:
    """Configuration for output files and directories."""
    output_dir: str = "test_results"
    figure_format: str = "pdf"
    figure_dpi: int = 300
    save_figures: bool = True
    show_figures: bool = False
    
    def __post_init__(self):
        """Create output directory if it doesn't exist."""
        if self.save_figures and not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print(f"Created output directory: {self.output_dir}")
    
    def get_path(self, filename: str) -> str:
        """Get full path for output file."""
        if not filename.endswith(f".{self.figure_format}"):
            filename = f"{filename}.{self.figure_format}"
        return os.path.join(self.output_dir, filename)
    
    def get_subdir(self, subdir: str) -> str:
        """Get or create subdirectory."""
        path = os.path.join(self.output_dir, subdir)
        if not os.path.exists(path):
            os.makedirs(path)
        return path


# Global configuration instance
_GLOBAL_CONFIG: Optional[OutputConfig] = None


def get_output_config() -> OutputConfig:
    """Get the global output configuration."""
    global _GLOBAL_CONFIG
    if _GLOBAL_CONFIG is None:
        _GLOBAL_CONFIG = OutputConfig()
    return _GLOBAL_CONFIG


def set_output_directory(path: str, figure_format: str = "pdf",
                        save_figures: bool = True, 
                        show_figures: bool = False) -> OutputConfig:
    """Set the global output configuration."""
    global _GLOBAL_CONFIG
    _GLOBAL_CONFIG = OutputConfig(
        output_dir=path,
        figure_format=figure_format,
        save_figures=save_figures,
        show_figures=show_figures
    )
    return _GLOBAL_CONFIG