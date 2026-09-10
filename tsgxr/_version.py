from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("tsgxr")
except PackageNotFoundError:
    # package is not installed
    pass
