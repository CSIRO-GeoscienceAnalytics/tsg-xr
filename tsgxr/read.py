import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytsg.parse_tsg
import xarray

# change a default setting for pytsg
# pytsg.parse_tsg.read_hires_dat = partial(pytsg.parse_tsg.read_hires_dat,  per_spectra=False)


def interpolate_section_depths(section_depths, ninterp):
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
        assert ninterp.dtype.kind in ["i"]
    return np.hstack(
        [
            np.linspace(mn, mx, nint)
            for (nint, (mn, mx)) in zip(
                ninterp,
                section_depths,
            )
        ]
    )


def cras_to_dataset(tsgdata, subsample=10):
    """
    Get the high resolution imagery from a TSG file, and optionally subsample it
    to a lower resolution.

    Parameters
    ----------
    tsgdata : pytsg.parse_tsg.TSG
        TSG dataset loaded with pytsg.

    Returns
    -------
    xarray.Dataset
        Dataset containing the RGB imagery.
    """
    section_depths = (
        tsgdata.nir.sampleheaders[["T", "L", "D"]]
        .apply(pd.to_numeric)
        .groupby(["T", "L"])
        .agg(["min", "max"])
        .values
    )

    depths = interpolate_section_depths(
        section_depths, [t.nlines for t in tsgdata.cras.section]
    )
    _dx = dy = np.median(np.diff(depths[:200]))
    horizontal = np.arange(0, tsgdata.cras.image.shape[1]) * dy
    horizontal -= horizontal.mean()
    cras = xarray.DataArray(
        tsgdata.cras.image,
        coords={
            "depth": ("sample", depths),
            "horizontal": horizontal,
            "channel": list("RGB"),
        },
        dims=("sample", "horizontal", "channel"),
        name="Image",
    )
    return cras[::subsample, ::subsample].to_dataset(name="RGB")


def load_tsg(
    directory,
    spectra="NIR",
    image=True,
    subsample_image=10,
    index_coord="sample",
    **kwargs,
):
    """
    Load a TSG dataset.

    Parameters
    ----------
    directory : str | pathlib.Path
        Directory of the TSG datset to load.
    spectra : str
        Which spectra to load by default, NIR or TIR.
    image : bool
        Whether to load the high-resolution RGB imagery.
    subsample_image : int
        Subsampling factor for the high-resolution RGB imagery.
        Default of 10 returns a 100x reduction in image size/
        1% of all pixels.
    index_coord : str
        Index coordinate to use for the dataset.
        Using "depth" requires some post-processing and dropping duplicates.

    Returns
    -------
    xarray.Dataset
        Dataset containing the spectra and associated data.
    """
    directory = Path(directory)
    tsgdata = pytsg.parse_tsg.read_package(directory, read_cras_file=image, **kwargs)
    dataset = tsg_to_xarray(tsgdata, spectra, index_coord=index_coord)
    if image:
        dataset["Image"] = cras_to_dataset(tsgdata, subsample=subsample_image)
        if (  # where evenly divisible by the subsample, we can assign sample-based coords
            124 % subsample_image == 0 and index_coord != "depth"
        ):
            pixels_per_sample = (
                dataset["Image"].coords["depth"].size / dataset.Spectra.sample.size
            )
            assert np.isclose(
                int(pixels_per_sample), pixels_per_sample
            )  # should be an even number
            pixels_per_sample = int(pixels_per_sample)
            # the image will have its own depth but otherwise the sample coords should transfer
            _crds = {
                k: ("sample", np.repeat(v.values, pixels_per_sample))
                if k != "sample"
                else np.repeat(v.values, pixels_per_sample)
                for k, v in dataset.Spectra.coords.items()
                if (k == "sample" or v.dims[0] == "sample") and k != "depth"
            }
            dataset["Image"]["RGB"] = dataset["Image"]["RGB"].assign_coords(_crds)
    return dataset


def tsg_to_xarray(tsgdata, spectra, index_coord="sample"):
    """
    Load a TSG spectral subset into Xarray.

    Parameters
    ----------
    tsgdata  : pytsg.parse_tsg.TSG
        TSG dataset loaded with pytsg.
    spectra : str
        Which spectra to load by default, NIR or TIR.
    index_coord : str
        Index coordinate to use for the dataset.
        Using "depth" requires some post-processing and dropping duplicates.

    Returns
    -------
    xarray.DataTree
        Dataset containing spectra and band headers.

    Todo
    -----
    * Consider indexing by depth instead of sample, after the fact.
    * Consider dropping Tray, Section, Depth (m) as they're duplicated as indexes.
    * Consider dropping SecDist (mm), TraySamp, SecSamp and NumFeats - they can be calculated.
    """
    spectraldata = getattr(tsgdata, spectra.lower())
    assert hasattr(spectraldata, "spectra"), (
        f"TSG Dataset does not have {spectra} data."
    )
    _coords = coords_from_sampleheaders(spectraldata)
    _sample_coords = {
        k: v
        for k, v in _coords.items()
        if k == "sample" or (isinstance(v, tuple) and v[0] == "sample")
    }

    scalar_data = spectraldata.scalars.copy()
    floatvals = scalar_data.select_dtypes(float).columns
    scalar_data[floatvals] = np.where(
        np.isclose(scalar_data.loc[:, floatvals].values, np.finfo("float32").min),
        np.nan,
        scalar_data.loc[:, floatvals.values],
    )
    # TODO: should spectra, products, Lidar be integrated into the same dataset?
    # could drop emtpy columns but is unlikely to be many
    products = scalar_data.set_index(
        pd.Series(scalar_data.index.values, name="sample")
    ).to_xarray()
    # TODO: some of these attributes could be propagated to other parts of the dataset
    products.attrs.update(
        {
            ch.name: [(i, v) for i, v in ch.classes.items()]
            for id, ch in spectraldata.classes.items()
        }
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
    #################################################################################
    # add the spectra, and move it to the top of the variable list
    spectra_ds = xarray.DataArray(
        spectraldata.spectra,
        coords=_coords,
        dims=("sample", "wavelength"),
    ).to_dataset(name="Spectra")
    #################################################################################
    # add the lidar data, sort out indexing
    if index_coord != "depth":
        # TODO: add depth as a secondary coordinate to this
        if tsgdata.lidar is not None:
            profilometer_ds = xarray.DataArray(
                tsgdata.lidar, coords={"sample": spectra_ds.sample.values}
            ).to_dataset(name="Lidar")
        else:
            profilometer_ds = None
        # alternate method for being able to index on depth for spectral without
        # dropping rows
        # specarr = specarr.set_xindex('depth')
    else:
        # remove samples where the depth is a duplicate, and sort by depth
        # to allow depth as an index
        fltr = pd.Series(spectra_ds.depth).duplicated().values
        spectra_ds = spectra_ds.sel(sample=~fltr)

        sortidx = np.argsort(spectra_ds.depth.values)
        # TODO: do we need to make sample a surrogate of depth here,
        # rather than depth being an independently indexed coordinate?
        spectra_ds["Spectra"] = spectra_ds["Spectra"][sortidx].swap_dims(
            {"sample": "depth"}
        )

        if tsgdata.lidar is not None:
            profilometer_ds = xarray.DataArray(
                tsgdata.lidar[~fltr][sortidx],
                coords={"depth": spectra_ds.depth.values},
            ).to_dataset(name="Lidar")
        else:
            profilometer_ds = None

    dataset = xarray.DataTree.from_dict(
        {"Spectra": spectra_ds, "Products": products.assign_coords(_sample_coords)}
    )
    if profilometer_ds is not None:
        dataset["Lidar"] = profilometer_ds.assign_coords(_sample_coords)
    return dataset


def coords_from_sampleheaders(spectraldata):
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
    sampleheaders = spectraldata.sampleheaders.rename(
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
        "wavelength": spectraldata.wavelength,
    }
    coords.update(
        {
            c: ("sample", d.values)
            for c, d in sampleheaders.items()
            if c not in ["sample"]
        }
    )
    return coords


def reorder_variables(
    ds,
    drop=[],  # ["Tray", "Section", "Depth (m)", "SecDist (mm)", "TraySamp", "SecSamp"],
    patterns=None,
):
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
