from pathlib import Path

import typer
import xarray
from tqdm.auto import tqdm

from . import __version__, find_TSG_datasets, load_tsg

app = typer.Typer()


@app.command()
def TSG2zarr(
    tsgdir: Path = typer.Argument(..., help="TSG directory."),
    output_dir: Path = typer.Option(
        None, help="Output directory for a Zarr data store."
    ),
    spectra: str = typer.Option("NIR", help="Whether to load NIR or TIR data."),
    index_coord: str = typer.Option(
        "Sample", help="Whether to use sample or depth as an index."
    ),
    image: bool = typer.Option(
        False, help="Whether to load the attached image, where available."
    ),
    subsample_image: int = typer.Option(
        10, help="Subsampling factor for imagery data."
    ),
):
    """
    Convert TSG file(s) to Zarr.
    """
    datasets = find_TSG_datasets(tsgdir)
    if datasets:
        for k, d in tqdm(datasets.items()):
            ds = load_tsg(
                d,
                spectra=spectra.upper(),
                image=image,
                subsample_image=subsample_image,
                index_coord=index_coord,
            )
            name = "_".join([str(ds.coords["hole"][0].values), spectra, ".zarr"])
            outdir = (
                (output_dir / name) if output_dir is not None else d / name
            )  # put it in the TSG directory if an output folder is not given
            ds.to_zarr(outdir, mode="w")  # overwrite if needed


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
    pass


if __name__ == "__main__":
    app()
