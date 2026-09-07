from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("tsg_xr")
except PackageNotFoundError:
    # package is not installed
    pass
