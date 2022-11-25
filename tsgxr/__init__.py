from . import _version

__version__ = _version.get_versions()["version"]

from .read import load_tsg
