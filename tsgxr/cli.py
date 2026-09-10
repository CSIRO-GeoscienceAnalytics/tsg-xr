import os
from pathlib import Path

os.environ["TYPER_USE_RICH"] = os.environ.get("TYPER_USE_RICH", "False")


import typer

try:
    from tqdm.auto import tqdm
except ImportError:

    def tqdm(iterable):
        return iterable


from . import __version__, find_TSG_datasets, open_tsg
from .util import Handle

logger = Handle(__name__, level="INFO")


app = typer.Typer()


@app.command()
def TSG2zarr(
    tsgdir: Path = typer.Argument(..., help="TSG directory or filepath."),  # noqa: B008
    output_dir: Path | None = typer.Option(  # noqa: B008
        None, help="Output directory for a Zarr data store."
    ),
    index_coord: str = typer.Option(
        "Sample", help="Whether to use sample or depth as an index."
    ),
    image: bool = typer.Option(
        False, help="Whether to load the attached image, where available."
    ),
    log: str = typer.Option(
        None,
        help="Whether to use logging, and if so what level (DEBUG, INFO)",
    ),
    zipfile: bool = typer.Option(
        True,
        "--zip",
        help="Whether to zip the Zarr archive upon creation.",
    ),
):
    """
    Convert TSG file(s) to Zarr.
    """
    if log is not None:
        Handle("tsgxr", level=log.upper())
    if (
        not tsgdir.is_dir() and tsgdir.suffix.lower() == ".tsg"
    ):  # pointing to a specific .tsg file
        logger.info(f"Loading TSG file: {tsgdir.name}")
        assert tsgdir.exists(), "Specified TSG file does not exist."
        datasets = {tsgdir.stem: tsgdir.parent}
    else:
        logger.info(f"Loading TSG files from directory: {tsgdir.name}")
        assert tsgdir.exists(), "Specified directory does not exist."
        datasets = find_TSG_datasets(tsgdir)

    if datasets:
        logger.info("Found datasets: {}".format(", ".join(k for k in datasets)))
    else:
        logger.warning(f"Found no datasets in {tsgdir.resolve()!s}.")

    if datasets:
        for k, d in tqdm(datasets.items()):
            ds = open_tsg(d, image=image, index_coord=index_coord)
            name = f"{ds.coords['hole'][0].values!s}{'.zip' if zipfile else '.zarr'}"
            # put it in the TSG directory if an output folder is not given
            outdir = output_dir if output_dir is not None else d
            logger.info(f"Creating Zarr archive {name} in {outdir!s}.")
            ds.to_zarr(outdir / name, mode="w")  # overwrite if needed


def version_callback(value: bool):
    if value:
        print(f"tsg-xr Version: {__version__}")
        raise typer.Exit()


@app.callback()
def callback(
    version: bool = typer.Option(
        None, "--version", callback=version_callback, is_eager=True
    ),
):
    """
    tsg-xr: A CLI tool for loading and transforming TSG files.
    """


if __name__ == "__main__":
    app()
