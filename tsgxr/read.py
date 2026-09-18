import re
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytsg.parse_tsg
import xarray

from .products import (
    get_available_systems,
    products_to_group_table,
    products_to_mineral_table,
)
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


def _subset_samples_no_depth_duplicates(
    da: xarray.DataArray | xarray.Dataset,
    template: xarray.Dataset | xarray.DataArray | None = None,
) -> xarray.DataArray:
    """
    Reindex a dataset's samples such that it's sorted by depth
    and depth-duplicates are removed.

    Parameters
    ----------
    da : xarray.DataArray | xarray.Dataset
        Dataset to reindex.
    template : xarray.DataArray | xarray.Dataset | None
        Template to use to get sample-depth indexing information from;
        where not provided this is expected to be taken from the input `da`.

    Returns
    -------
    xarray.DataArray
        Index of samples fitting the requirements.

    Notes
    -----
    This is a lossy process where depth is duplicated.
    """
    if template is None:
        template = da
    fltr = pd.Series(template.depth).duplicated().values
    return (
        da.sample.sel(sample=~fltr)
        .isel(sample=np.argsort(template.sel(sample=~fltr).depth.values))
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
    arrangement = []
    # these occur where collapse_products=True, if they're there put them first.
    arrangement += sorted([v for v in ds.data_vars if "_Grp" in v or "_Min" in v])
    # The first of these are e.g. coordinates and would typically get dropped
    arrangement += [
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
        if v in ds.data_vars and v not in arrangement
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
        drop_vars = [
            "HoleID",
            "Date",
            "Depth (m)",
            "Tray",
            "Section",
            "NumFeats",
            "SecDist (mm)",
            "SecSamp",
            "TraySamp",
        ]
    scalar_data = scalars.copy()  # pd.DataFrame
    if drop_vars:
        scalar_data = scalar_data.drop(
            columns=[v for v in drop_vars if v in scalar_data.columns]
        )
        # some of the dropped variables are classes, which might also have colormaps
        scalar_data.attrs = {
            k: v
            for k, v in scalar_data.attrs.items()
            if ((k not in drop_vars) and (k.replace("_Colors", "") not in drop_vars))
        }
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
        # # [f.replace(grp, "") for f in arr["feature"].values]
        # NOTE: these are just a sequence
        arr["feature"] = np.arange(arr["feature"].size, dtype=np.uint8)
        # NOTE: numbers of features are given in the metadata, so could
        # use that rather than check against zero here
        arr = xarray.where(arr == 0, np.nan, arr)
        arr.attrs = {}

        products[grp + "s"] = arr.T  # put depth/sample along the first axis
    if drop_singular:
        # drop variables which only have one value; these are typically 0, 1, nan or 'Default
        products = products.drop_vars(
            [
                k
                for k in products.data_vars
                if products[k].dims == ("sample",)
                and pd.unique(pd.Series(products[k])).size == 1
            ]
        )
    if collapse_products:
        systems = get_available_systems(products)
        dropprod = []
        for system in systems:
            # NOTE: this will also drop error, NIL_Stat, SNR
            dropprod += [k for k in products.data_vars if k.endswith(system)]
            products[f"{system}_Grp"] = (
                products_to_group_table(products, which=system)
                .to_xarray()
                .to_dataarray(f"{system}group")
                .T
            )
            products[f"{system}_Min"] = (
                products_to_mineral_table(products, which=system)
                .to_xarray()
                .to_dataarray(f"{system}mineral")
                .T
            )
        products = products.drop_vars(dropprod)
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
    collapse_products: bool = False,
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
    collapse_products : bool
        Whether to collapse products to a singular table per system,
        rather than multiple e.g. Min1, Min2, ..

    Returns
    -------
    xarray.DataTree
        DataTree containing the spectra and associated data.
    """
    directory = Path(directory)

    lidar = next(directory.glob("*tsg_hires.dat*"))

    files = sorted(
        [f for f in directory.glob("*.bip*") if "cras" not in f.stem],
        key=lambda x: x.stem,
    ) + ([lidar] if lidar else [])

    def _load(fpath):
        if fpath.suffix == ".bip":
            ds = xarray.open_dataset(
                fpath,
                engine="tsg",
                collapse_products=collapse_products,
            ).drop_vars("half")
            if not lazy:
                ds = ds.load()
        else:
            # TODO: lazy loader? Probably not worth it here
            ds = xarray.DataArray(
                pytsg.parse_tsg.read_hires_dat(lidar), dims=("sample",)
            ).to_dataset(name="Lidar")
        return ds

    D = {
        SPECTRAL_MAPPING.get(f.stem.split("_")[-1])
        if f.suffix == ".bip"
        else "Lidar": v
        for f, v in zip(
            files,
            joblib.Parallel(
                backend="threading", n_jobs=len(files), return_as="generator"
            )(joblib.delayed(_load)(fpath) for fpath in files),
        )
    }
    template = D[next(iter([k for k in D if "Spectra" in D[k].data_vars]))]
    if index_coord == "depth":
        # TODO: can we reindex depth once, for all of these?
        samples = _subset_samples_no_depth_duplicates(template)
        for k in D:
            if k == "Lidar":  # Lidar still needs the extra coordinates
                D[k] = D[k].assign_coords(
                    {
                        k: v
                        for k, v in template.coords.items()
                        if k == "sample" or (isinstance(v, tuple) and v[0] == "sample")
                    }
                )
            D[k] = D[k].sel(sample=samples).swap_dims({"sample": "depth"})
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
            # NOTE: these uses the template spectral dataset accessed above
            # this version hasn't yet been depth-reindexed!
            depths = interpolate_section_depths(
                template, image_ds.Image.attrs["section"]
            )
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
                    image_ds.coords["depth"].size
                    / template["Spectra"].coords["sample"].size
                )
                assert np.isclose(
                    int(pixels_per_sample), pixels_per_sample
                )  # should be an even number
                pixels_per_sample = int(pixels_per_sample)
                # the image will have its own depth but otherwise the sample coords should transfer
                image_ds = image_ds.assign_coords(
                    {
                        k: ("x", np.repeat(v.values, pixels_per_sample))
                        for k, v in template["Spectra"].coords.items()
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
    del template
    return xarray.DataTree.from_dict(D)
