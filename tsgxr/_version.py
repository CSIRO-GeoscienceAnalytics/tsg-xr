from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("tsg_xr")
except PackageNotFoundError:
    # package is not installed
    pass
