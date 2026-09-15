import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytsg.parse_tsg
import xarray

from .util import Handle

logger = Handle(__name__)

SPECTRAL_MAPPING = {"tsg": "NIR", "tir": "TIR", "mir": "MIR"}


def interpolate_section_depths(
    template: xarray.DataArray | xarray.Dataset, sections: np.rec.recarray
) -> np.ndarray:
    """
    Interpolate section depths based on a number of interpolated samples.

    Parameters
    ----------
    template: xarray.DataArray | xarray.Dataset
        Data structure containing the requisite tray, section and depth coordinates.
    sections : np.rec.recarray
        Record array containing the 'nlines' information per-section.

    Returns
    -------
    numpy.ndarray (nsamples, )
        Interpolated within-section depths.
    """

    idx = pd.MultiIndex.from_arrays(
        [template.tray.values, template.section.values], names=["tray", "section"]
    )
    section_depths = (
        pd.Series(template.depth.values).groupby(idx).agg(["min", "max"]).to_xarray()
    )
    section_depths["index"] = pd.MultiIndex.from_arrays(
        np.array(list(section_depths["index"].values)).T, names=["tray", "section"]
    )

    assert len(sections) == len(section_depths["min"]), (
        f"Length of sections record array ({len(sections)}) doesn't match the depth metadata ({len(section_depths['min'])})."
    )

    if (sections["nlines"] == sections["nlines"][0]).all():
        # all have the same number of lines
        section_depth_interp = (
            np.linspace(0, 1, sections["nlines"][0], dtype=section_depths["max"].dtype)[
                None, :
            ]
            * (section_depths["max"] - section_depths["min"]).values[:, None]
            + section_depths["min"].values[:, None]
        ).ravel()
    else:
        section_depth_interp = np.hstack(
            [
                np.linspace(mn, mx, nint, dtype=section_depths["max"].dtype)
                for (mn, mx, nint) in zip(
                    section_depths["min"].values,
                    section_depths["max"].values,
                    sections["nlines"],
                    strict=True,
                )
            ]
        )
    return section_depth_interp


def _reindex_depth(
    da: xarray.DataArray | xarray.Dataset,
    template: xarray.Dataset | xarray.DataArray | None = None,
) -> xarray.DataArray | xarray.Dataset:
    """
    Reindex a dataset such that it's indexed by depth.

    Parameters
    ----------
    da : xarray.DataArray | xarray.Dataset
        Dataset to reindex.
    template : xarray.DataArray | xarray.Dataset | None
        Template to use to get sample-depth indexing information from;
        where not provided this is expected to be taken from the input `da`.

    Returns
    -------
    xarray.DataArray | xarray.Dataset
        Reindexed dataset.

    Notes
    -----
    This is a lossy process where depth is duplicated.
    """
    # remove samples where the depth is a duplicate, and sort by depth
    # to allow depth as an index
    if template is None:
        template = da
    fltr = pd.Series(template.depth).duplicated().values
    return (
        da.sel(sample=~fltr)
        .isel(sample=np.argsort(template.sel(sample=~fltr).depth.values))
        .swap_dims({"sample": "depth"})
        .sortby("depth")
    )


def reorder_variables(
    ds: xarray.Dataset,
    drop: list
    | None = None,  # ["Tray", "Section", "Depth (m)", "SecDist (mm)", "TraySamp", "SecSamp"],
    patterns: list | None = None,
) -> xarray.Dataset:
    """
    Reorder the variables within an Xarray dataset containing TSG data such that
    it's more easily visually navigated (note this does not persist upon serialization).

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to reorder.
    drop : list
        Variables to exclude (typically either duplciated in indexes or easily calculated).
    patterns : list
        List of regex patterns to match column groups of band headers within TSG files.

    Returns
    -------
    ds : xarray.Dataset
        Reordered dataset.
    """
    if drop is None:
        drop = []
    if patterns is None:
        patterns = [
            r"Grp\d*",
            r"Min\d*",
            r"Wt\d*",
            r"Error\d*",
            "SNR",
            "NIL_Stat",
            "Cust",
            "Bound_Water",
            "Unbound_Water",
        ]
    arrangement = [
        v
        for v in [
            "HoleID",
            "Date",
            "Depth (m)",
            "Tray",
            "Section",
            "Centres",
            "Depths",
            "Widths",
        ]
        if v in ds.data_vars
    ]

    arrangement += [
        v
        for ptn in patterns
        for v in sorted([v for v in ds.data_vars if re.match(ptn, v)])
    ]

    others = sorted(
        set(
            [v for v in ds.data_vars if (v not in arrangement)]
            + ["Flags"]
            + [v for v in ds.data_vars if v.lower() == v]
        )
    )
    ds = ds[
        [v for v in arrangement + others if (v not in drop) and (v in ds.data_vars)]
    ]
    return ds


def product_dataset_to_xarray(
    scalars: pd.DataFrame,
    collapse_products: bool = False,
    drop_vars=None,
    drop_singular=True,
) -> xarray.Dataset:
    """
    Transform a set of spectral products/scalars into xarray.align

    Parameters
    ----------
    scalars  : pandas.DataFrame
        Dataframe of loaded by pytsg.
    collapse_products : bool
        Whether to collapse products to a singular table per system,
        rather than multiple e.g. Min1, Min2, ..
    drop_vars : list
        List of variables to drop, typically due to being duplicated as coordinates.
    drop_singular : bool
        Whether to drop singular values which are propagated across the products.

    Returns
    -------
    xarray.Dataset

    Notes
    -----
    * Note that group is essentially redundant, could be a coordinate on mineral.
    """
    if drop_vars is None:
        drop_vars = ["HoleID", "Date", "Depth (m)", "Tray", "Section"]
    scalar_data = scalars.copy()  # pd.DataFrame
    if drop_vars:
        scalar_data = scalar_data.drop(
            columns=[v for v in drop_vars if v in scalar_data.columns]
        )
    # could drop emtpy columns but is unlikely to be many
    scalar_data.index.name = "sample"
    products = scalar_data.to_xarray()
    products.attrs.update(scalar_data.attrs)  # propagate attributes
    for grp in ["Centre", "Depth", "Width"]:
        arr = (
            products[[v for v in products.data_vars if re.match(grp + r"\d+", v)]]
            .to_array()
            .rename({"variable": "feature"})
        )
        products = products[
            [v for v in products.data_vars if v not in arr.coords["feature"]]
        ]
        arr["feature"] = [f.replace(grp, "") for f in arr["feature"].values]
        arr = xarray.where(arr == 0, np.nan, arr)
        arr.attrs = {}

        products[grp + "s"] = arr
    if drop_singular:
        # drop variables which only have one value; these are typically 0, 1, nan or 'Default
        products = products.drop_vars(
            [
                k
                for k in products.data_vars
                if products[k].dims == ("depth",)
                and pd.unique(pd.Series(products[k])).size == 1
            ]
        )
    # convert traynames, otherwise occasionally converted to integers
    if "Tray" in products:
        products["Tray"] = products["Tray"].astype("<U16")
    products = reorder_variables(products)
    return products


def open_tsg(
    directory: str | Path,
    image: bool = True,
    index_coord: str = "sample",
    lazy: bool = True,
    chunks: int | dict | None = None,
    **kwargs,
) -> xarray.DataTree:
    """
    Open a TSG dataset.

    Parameters
    ----------
    directory : str | pathlib.Path
        Directory of the TSG datset to load.
    image : bool
        Whether to load the high-resolution RGB imagery.
    index_coord : str
        Index coordinate to use for the dataset.
        Using "depth" requires some post-processing and dropping duplicates.
    lazy : bool
        Whether to load the dataset lazily (default), or otherwise
        load after the data structure is ready.
    chunks : int | dict | None
        Chunking specification for the dataset, if you're planning to
        use `dask`.

    Returns
    -------
    xarray.DataTree
        DataTree containing the spectra and associated data.
    """
    directory = Path(directory)

    # TODO: read each of the file pairs using the xarray backend

    spectral_bips = [f for f in directory.glob("*.bip*") if "cras" not in f.stem]
    D = {}
    for f in spectral_bips:
        sset = SPECTRAL_MAPPING.get(f.stem.split("_")[-1])
        ds = xarray.open_dataset(f, engine="tsg").drop_vars("half")
        if not lazy:
            ds = ds.load()
        D = {**D, f"{sset}": ds}
    lidar = next(directory.glob("*tsg_hires.dat*"))
    if lidar:
        prof_da = xarray.DataArray(  # TODO: lazy loader?
            pytsg.parse_tsg.read_hires_dat(lidar), dims=("sample",)
        ).assign_coords(
            {
                k: v
                for k, v in ds.coords.items()
                if k == "sample" or (isinstance(v, tuple) and v[0] == "sample")
            }
        )
        if index_coord == "depth":
            prof_da = _reindex_depth(prof_da, template=ds)
        if chunks:
            prof_da = prof_da.chunk(
                chunks
                if isinstance(chunks, int)
                else {k: v for k, v in chunks.items() if k in prof_da.dims}
            )
        D["Lidar"] = prof_da.to_dataset(name="Lidar")

    if index_coord == "depth":  # TODO: rechunk?
        for k in D:
            if "Spectra" in D[k].data_vars:
                D[k] = _reindex_depth(D[k])
                if chunks:
                    D[k] = D[k].chunk(
                        chunks
                        if isinstance(chunks, int)
                        else {k: v for k, v in chunks.items() if k in D[k].dims}
                    )

    if image:
        crasfile = list(directory.glob("*cras.bip*"))
        if crasfile:
            crasfile = crasfile[0]
        if not crasfile:
            raise FileNotFoundError(
                "CRAS file matching *cras.bip not found in directory."
            )
        else:
            image_ds = xarray.open_dataset(
                crasfile, engine="lazycras" if lazy else "cras"
            )
            # NOTE: these uses the last spectral dataset accessed above
            # this version hasn't yet been depth-reindexed!
            depths = interpolate_section_depths(ds, image_ds.Image.attrs["section"])
            # TODO: assign section-part IDs?
            _dx = dy = np.median(np.diff(depths[:200]))
            horizontal = np.arange(0, image_ds.Image.shape[1]) * dy
            horizontal -= horizontal.mean()
            image_ds = image_ds.assign_coords(
                depth=("x", depths), width=("y", horizontal)
            )

            if index_coord == "depth":
                # NOTE: width only appears here to make aspect-equal plotting easier - it's either pixel-pixel or depth-width
                # there are no duplicate depths in the image, so we dont' need to deduplicate this
                # TODO: assign sample-based coordinates as per depth ranges in the deduplicated sample spectra image
                image_ds = (
                    image_ds.swap_dims({"x": "depth"})
                    .swap_dims({"y": "width"})
                    .sortby("depth")
                    .assign_coords(
                        {
                            "section": ("depth", image_ds.section.values),
                            "tray": ("depth", image_ds.tray.values),
                        }
                    )
                )

            else:
                #  we can assign sample-based coords
                pixels_per_sample = (
                    image_ds.coords["depth"].size / ds["Spectra"].coords["sample"].size
                )
                assert np.isclose(
                    int(pixels_per_sample), pixels_per_sample
                )  # should be an even number
                pixels_per_sample = int(pixels_per_sample)
                # the image will have its own depth but otherwise the sample coords should transfer

                image_ds = image_ds.assign_coords(
                    {
                        k: ("x", np.repeat(v.values, pixels_per_sample))
                        for k, v in ds["Spectra"].coords.items()
                        if (k == "sample" or v.dims[0] == "sample") and k != "depth"
                    }
                )

            if chunks:
                image_ds = image_ds.chunk(
                    chunks
                    if isinstance(chunks, int)
                    else {k: v for k, v in chunks.items() if k in image_ds.dims}
                )
            D["Image"] = image_ds
    return xarray.DataTree.from_dict(D)
