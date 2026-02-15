from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("hyspxr")
except PackageNotFoundError:
    # package is not installed
    pass
