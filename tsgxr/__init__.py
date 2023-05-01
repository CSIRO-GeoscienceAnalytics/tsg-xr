from pathlib import Path 

from . import _version

__version__ = _version.get_versions()["version"]

from .read import load_tsg


def find_TSG_datasets(parent_directory):
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
            set(list(Path(parent_directory).glob("**/*_tsg.tsg")) + list(Path(parent_directory).glob("**/*_tsg_tir.tsg")))
        )
    }
