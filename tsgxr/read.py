import inspect
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
    section_depths: np.ndarray, ninterp: int | np.ndarray
) -> np.ndarray:
    """
    Interpolate section depths based on a number of interpolated samples.

    Parameters
    ----------
    section_depths : numpy.ndarray (n_sections, 2)
        Array containing the minimum and maximum depths of each section.
    ninterp : int | numpy.ndarray (n_sections)
        Either a constant number of divisions or a number of divisions per section.

    Returns
    -------
    numpy.ndarray (nsamples, )
        Interpoalted within-section depths.
    """

    assert isinstance(ninterp, (int, np.ndarray, list))
    if isinstance(ninterp, int):
        ninterp = np.ones(section_depths.shape[0]) * ninterp
    else:
        ninterp = np.array(ninterp)
        assert ninterp.dtype.kind in ["i", "u"]
    return np.hstack(
        [np.linspace(mn, mx, nint) for (nint, (mn, mx)) in zip(ninterp, section_depths)]
    )


def coords_from_sampleheaders(headers: pd.DataFrame, wavelengths) -> dict:
    """
    Turn the sample headers of a TSG spectral subset into coordinates.

    Parameters
    ----------
    spectraldata  : pytsg.parse_tsg.Spectra
        Spectral subset loaded with pytsg.

    Returns
    -------
    coords : dict
        Mapping of coordinate names to values, and in the case of non-index coordinates
        the corresponding index coordinate.
    """
    sampleheaders = (headers).rename(
        columns={
            "sample": "sample",
            "T": "tray",
            "L": "section",
            "P": "section-part",
            "D": "depth",
            "X": "section-position",
            "H": "hole",
        }
    )
    for k in sampleheaders.columns:  # try to convert numeric data
        try:
            sampleheaders[k] = sampleheaders[k].apply(pd.to_numeric)
        except ValueError:
            pass

    # note that depths can be duplicated, so would need to be
    # post-processed to be used as an index
    coords = {
        "sample": sampleheaders["sample"].values,
        "wavelength": ("band", wavelengths),
    }
    coords.update(
        {
            c: ("sample", d.values)
            for c, d in sampleheaders.items()
            if c not in ["sample"]
        }
    )
    return coords


def _reindex_depth(
    da: xarray.DataArray, template: xarray.Dataset | xarray.DataArray | None = None
) -> xarray.DataArray:

    # remove samples where the depth is a duplicate, and sort by depth
    # to allow depth as an index
    if template is None:
        template = da
    fltr = pd.Series(template.depth).duplicated().values
    template = template.sel(sample=~fltr)
    sortidx = np.argsort(template.depth.values)
    return da.isel(sample=sortidx).swap_dims({"sample": "depth"}).sortby("depth")


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


def product_dataset_to_xarray(scalars: pd.DataFrame, classes: dict) -> xarray.Dataset:
    """
    Transform a set of spectral products/scalars into xarray.align

    Parameters
    ----------
    scalars  : pandas.DataFrame
        Dataframe of loaded by pytsg.
    classes : dict
        Mapping of classes, to be added as attributes.

    Returns
    -------
    xarray.Dataset
    """
    scalar_data = scalars.copy()  # pd.DataFrame
    floatvals = scalar_data.select_dtypes(float).columns
    scalar_data[floatvals] = np.where(
        np.isclose(scalar_data.loc[:, floatvals].values, np.finfo("float32").min),
        np.nan,
        scalar_data.loc[:, floatvals.values],
    )
    # could drop emtpy columns but is unlikely to be many
    products = scalar_data.set_index(
        pd.Series(scalar_data.index.values, name="sample")
    ).to_xarray()
    products.attrs.update(
        {ch.name: [(i, v) for i, v in ch.classes.items()] for ID, ch in classes.items()}
    )
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
    # convert traynames, otherwise occasionally converted to integers
    products["Tray"] = products["Tray"].astype("<U16")
    products = reorder_variables(products)
    return products


def spectral_dataset_to_xarray(
    spectra: pytsg.parse_tsg.Spectra,
    index_coord: str = "sample",
    chunks: dict | int | None = None,
) -> xarray.Dataset:
    """
    Load a TSG spectral subset into Xarray.

    Parameters
    ----------
    spectra  : pytsg.parse_tsg.Spectra
        TSG spectral dataset loaded with pytsg.
    index_coord : str
        Index coordinate to use for the dataset.
        Using "depth" requires some post-processing and dropping duplicates.
    chunks : int | dict
        Chunking to use for the dataset.

    Returns
    -------
    xarray.Dataset
    """
    _coords = coords_from_sampleheaders(spectra.sample_headers, spectra.wavelength)
    _sample_coords = {
        k: v
        for k, v in _coords.items()
        if k == "sample" or (isinstance(v, tuple) and v[0] == "sample")
    }
    # handle products
    products = product_dataset_to_xarray(spectra.scalars, spectra.classes)
    ###############################################################################s##
    # add the spectra, and move it to the top of the variable list
    spectra_da = xarray.DataArray(
        spectra.spectra,
        coords=_coords,
        dims=("sample", "band"),
    )
    if index_coord == "depth":
        # have to do products first, otherwise spectra_da is changed
        products = _reindex_depth(products, template=spectra_da)
        spectra_da = _reindex_depth(spectra_da)

    if chunks:
        spectra_da = spectra_da.chunk(
            chunks
            if isinstance(chunks, int)
            else {k: v for k, v in chunks.items() if k in spectra_da.dims}
        )
        products = products.chunk(
            chunks
            if isinstance(chunks, int)
            else {k: v for k, v in chunks.items() if k in products.dims}
        )

    return products.assign_coords(_sample_coords).assign(Spectra=spectra_da)


def tsg_to_xarray(
    tsgdata,
    index_coord: str = "sample",
    chunks: dict | int | None = None,
) -> xarray.DataTree:
    """
    Load an entire TSG dataset into Xarray.

    Parameters
    ----------
    tsgdata  : pytsg.parse_tsg.TSG
        TSG dataset loaded with pytsg.
    index_coord : str
        Index coordinate to use for the dataset.
        Using "depth" requires some post-processing and dropping duplicates.
    chunks : int | dict
        Chunking to use for the dataset.

    Returns
    -------
    xarray.DataTree
        Data tree containing spectra and band headers.

    Todo
    -----
    * Consider dropping Tray, Section, Depth (m) as they're duplicated as indexes.
    * Consider dropping SecDist (mm), TraySamp, SecSamp and NumFeats - they can be calculated.
    """
    if chunks is None:
        chunks = {}

    DT = xarray.DataTree()

    for subset in ["nir", "mir", "tir"]:
        if hasattr(tsgdata, subset) and not inspect.isclass(getattr(tsgdata, subset)):
            specds: xarray.Dataset = spectral_dataset_to_xarray(
                getattr(tsgdata, subset), index_coord=index_coord, chunks=chunks
            )
            DT[subset.upper()] = xarray.DataTree.from_dict(
                {
                    "Spectra": specds["Spectra"].to_dataset(name="Spectra"),
                    "Products": specds.drop_vars("Spectra").drop_dims("band"),
                }
            )

    #################################################################################
    # add the lidar data, sort out indexing
    if tsgdata.lidar is not None:
        # the coordinates used here need to be the sample ones
        prof_da = xarray.DataArray(tsgdata.lidar, dims=("sample",)).assign_coords(
            {
                k: v
                for k, v in coords_from_sampleheaders(
                    getattr(tsgdata, subset).sample_headers,
                    getattr(tsgdata, subset).wavelength,
                ).items()
                if k == "sample" or (isinstance(v, tuple) and v[0] == "sample")
            }
        )
        if index_coord == "depth":
            prof_da = _reindex_depth(prof_da)
        if chunks:
            prof_da = prof_da.chunk(
                chunks
                if isinstance(chunks, int)
                else {k: v for k, v in chunks.items() if k in prof_da.dims}
            )
        DT["Lidar"] = prof_da.chunk().to_dataset(name="Lidar")

    return DT


def open_tsg(
    directory,
    image=True,
    index_coord="sample",
    lazy=True,
    chunks=None,
    **kwargs,
):
    """
    Open a TSG dataset.

    Parameters
    ----------
    directory : str | pathlib.Path
        Directory of the TSG datset to load.
    spectra : str
        Which spectra to load by default, NIR or TIR.
    image : bool
        Whether to load the high-resolution RGB imagery.
    index_coord : str
        Index coordinate to use for the dataset.
        Using "depth" requires some post-processing and dropping duplicates.

    Returns
    -------
    xarray.Dataset
        Dataset containing the spectra and associated data.
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
        D = {
            **D,
            f"{sset}/Spectra": ds[["Spectra"]],
            f"{sset}/Products": ds.drop_vars("Spectra"),
        }
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
            if "Products" in k:
                D[k] = _reindex_depth(D[k], template=D[f"{k.split('/')[0]}/Spectra"])
                if chunks:
                    D[k] = D[k].chunk(
                        chunks
                        if isinstance(chunks, int)
                        else {k: v for k, v in chunks.items() if k in D[k].dims}
                    )
        for k in D:
            if "Spectra" in k:
                D[k] = _reindex_depth(D[k])
                if chunks:
                    D[k] = D[k].chunk(
                        chunks
                        if isinstance(chunks, int)
                        else {k: v for k, v in chunks.items() if k in D[k].dims}
                    )

    DT = xarray.DataTree.from_dict(D)
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
            section_depths = (
                (
                    ds.depth.groupby(["tray", "section"]).map(
                        lambda x: xarray.Dataset({"min": x.min(), "max": x.max()})
                    )
                )
                .stack(s=("tray", "section"))
                .to_dataarray("metric")
                .T.dropna(how="any", dim="s")
                .values
            )
            depths = interpolate_section_depths(
                section_depths, [t.nlines for t in image_ds.Image.attrs["section"]]
            )
            _dx = dy = np.median(np.diff(depths[:200]))
            horizontal = np.arange(0, image_ds.Image.shape[1]) * dy
            horizontal -= horizontal.mean()
            image_ds = image_ds.assign_coords(
                depth=("x", depths), width=("y", horizontal)
            )

            if index_coord == "depth":
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
            DT["Image"] = image_ds
    return DT
