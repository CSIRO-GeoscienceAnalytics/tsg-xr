from pathlib import Path

from ._version import __version__
from .read import open_tsg
from .util import Handle

logger = Handle(__name__)


__all__ = ["__version__", "find_TSG_datasets", "open_tsg"]


def find_TSG_datasets(parent_directory: Path | str) -> dict:
    """
    Check a directory for subdirectories containing Hylogger TSG datasets.

    Parameters
    ----------
    parent_directory : str | pathlib.Path
        Directory containing Hylogger datasets.

    Returns
    -------
    dict
        Dictionary mapping of hole names to subdirectories containing
        Hylogger datasets (the directory to be passed to `pytsg` or
        `tsg-xr`).
    """
    return {
        fpath.stem.replace("_tsg", ""): fpath.parent
        for fpath in sorted(
            set(
                list(Path(parent_directory).glob("**/*_tsg.tsg"))
                + list(Path(parent_directory).glob("**/*_tsg_tir.tsg"))
            )
        )
    }
