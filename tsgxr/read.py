from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import xarray

import pytsg.parse_tsg


def load_lidar(hirespath, per_spectra=True):
    """
    Patch for loading high resolution profilometer data
    with pytsg.

    Parameters
    ----------
    hirespath : str | pathlib.path
        Path to the hires.dat file containing profilometer data.
    per_spectra : bool
        Whether to aggregate the profilometer data to a per-sample
        /per-spectra basis, allowing direct integration.

    Returns
    -------
    numpy.ndarray
        Profilometer data array.
    """
    with open(hirespath, "rb") as f:
        idchar = f.read(20)  # CoreLog high-res 1.0
        nsclr, nl, nsps = np.fromfile(f, np.int32, 3)
        minp, maxp = np.fromfile(f, np.float32, 2)
        _prof = f.read(len("Profilometer"))
        # four int8/bytes, two int16 or one int32 gives 64 * 8 * n_spectra profilometer samples
        _ = np.fromfile(f, np.ubyte, 4)
        assert (_ == 0).all()
        lidar = np.fromfile(f, dtype=np.float32)
        # could do additional validation here
        lidar[lidar < minp] = np.nan

    if per_spectra:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            lidar = np.nanmean(lidar.reshape(-1, nsps), axis=1)
    return lidar


pytsg.parse_tsg.read_hires_dat = load_lidar  # patch the function


def load_tsg(
    directory, spectra="NIR", image=True, subsample_image=10, index_coord="sample"
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
        Dataset containing the spectra and assocaited data.
    """
    directory = Path(directory)
    tsgdata = pytsg.parse_tsg.read_package(directory, read_cras_file=image)
    dataset = tsg_to_xarray(tsgdata, spectra, index_coord=index_coord)
    if image:
        dataset["Image"] = cras_to_dataarray(tsgdata, subsample=subsample_image)
    dataset = reorder_variables(dataset)
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
    xarray.Dataset
        Dataset containing spectra and band headers.

    Todo
    -----
    * Consider indexing by depth instead of sample, after the fact.
    * Consider dropping Tray, Section, Depth (m) as they're duplicated as indexes.
    * Consider dropping SecDist (mm), TraySamp, SecSamp and NumFeats - they can be calculated.
    """
    spectraldata = getattr(tsgdata, spectra.lower())
    scalar_data = spectraldata.scalars.copy()

    scalar_data = spectraldata.scalars.copy()
    floatvals = scalar_data.select_dtypes(float).columns
    scalar_data[floatvals] = np.where(
        np.isclose(scalar_data.loc[:, floatvals].values, np.finfo("float32").min),
        np.nan,
        scalar_data.loc[:, floatvals.values],
    )
    # could drop emtpy columns but is unlikely to be many
    dataset = scalar_data.set_index(
        pd.Series(scalar_data.index.values, name="sample")
    ).to_xarray()
    dataset.attrs.update(
        {
            ch.name: [(i, v) for i, v in ch.classes.items()]
            for id, ch in spectraldata.classes.items()
        }
    )

    for grp in ["Centre", "Depth", "Width"]:
        arr = (
            dataset[[v for v in dataset.data_vars if re.match(grp + "\d+", v)]]
            .to_array()
            .rename({"variable": "feature"})
        )
        dataset = dataset[
            [v for v in dataset.data_vars if v not in arr.coords["feature"]]
        ]
        arr["feature"] = [f.replace(grp, "") for f in arr["feature"].values]
        arr = xarray.where(arr == 0, np.nan, arr)
        arr.attrs = {}

        dataset[grp + "s"] = arr

    # convert traynames, otherwise occasionally converted to integers
    dataset["Tray"] = dataset["Tray"].astype("<U16")
    # add the spectra, and move it to the top of the variable list
    coords = coords_from_sampleheaders(spectraldata)
    specarr = xarray.DataArray(
        spectraldata.spectra, coords=coords, dims=("sample", "wavelength")
    )
    if index_coord != "depth":
        profilometer = xarray.DataArray(
            tsgdata.lidar, coords={"sample": specarr.sample.values}
        )
        # alternate method for being able to index on depth for spectral without
        # dropping rows
        # specarr = specarr.set_xindex('holedepth')
    else:
        # remove samples where the depth is a duplicate, and sort by depth
        # to allow depth as an index
        fltr = pd.Series(specarr.holedepth).duplicated().values
        specarr = specarr.sel(sample=~fltr)

        sortidx = np.argsort(specarr.holedepth.values)
        specarr = specarr[sortidx].swap_dims({"sample": "holedepth"})

        profilometer = xarray.DataArray(
            tsgdata.lidar[~fltr][sortidx],
            coords={"holedepth": specarr.holedepth.values},
        )
    dataset["Spectra"] = specarr
    dataset["Lidar"] = profilometer
    dataset = dataset[["Spectra"] + [v for v in dataset.data_vars if v != "Spectra"]]
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
    sampleheaders = spectraldata.sampleheaders.apply(
        pd.to_numeric, errors="ignore"
    ).rename(
        columns={
            "sample": "sample",
            "T": "tray",
            "L": "section",
            "P": "section-part",
            "D": "holedepth",  # holedepth here as we can't easily deal with two depth indexes
            "X": "section-position",
            "H": "hole",
        }
    )
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


def cras_to_dataarray(tsgdata, subsample=10):
    """
    Get the high resolution imagery from a TSG file, and optionally subsample it
    to a lower resolution.

    Parameters
    ----------
    tsgdata : pytsg.parse_tsg.TSG
        TSG dataset loaded with pytsg.

    Returns
    -------
    xarray.DataArray
        Array containing the RGB imagery.
    """
    depths = np.hstack(
        [
            np.linspace(mn, mx, t.nlines)
            for (t, (mn, mx)) in zip(
                tsgdata.cras.section,
                tsgdata.nir.sampleheaders[["T", "L", "D"]]
                .apply(pd.to_numeric)
                .groupby(["T", "L"])
                .agg(["min", "max"])
                .values,
            )
        ]
    )
    dx = dy = np.median(np.diff(depths[:200]))
    horizontal = np.arange(0, tsgdata.cras.image.shape[1]) * dy
    horizontal -= horizontal.mean()
    cras = xarray.DataArray(
        tsgdata.cras.image,
        coords={"depth": depths, "horizontal": horizontal, "channel": list("RGB")},
    )
    return cras[::subsample, ::subsample]


def reorder_variables(
    ds,
    drop=["Tray", "Section", "Depth (m)", "SecDist (mm)", "TraySamp", "SecSamp"],
    patterns=[
        "Grp\d*",
        "Min\d*",
        "Wt\d*",
        "Error\d*",
        "SNR",
        "NIL_Stat",
        "Cust",
        "Bound_Water",
        "Unbound_Water",
    ],
):
    """
    Reorder the variables within an Xarray dataset containing TSG data such that
    it's more easily visually navigated (note this does not persist upon
    serialization).

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
    arrangement = [
        "HoleID",
        "Date",
        "Depth (m)",
        "Tray",
        "Section",
        "Spectra",
        "Image",
        "Centres",
        "Depths",
        "Widths",
    ]

    arrangement += [
        v
        for ptn in patterns
        for v in sorted([v for v in ds.data_vars if re.match(ptn, v)])
    ]

    others = sorted(
        list(
            set(
                [v for v in ds.data_vars if (v not in arrangement)]
                + ["Flags"]
                + [v for v in ds.data_vars if v.lower() == v]
            )
        )
    )
    ds = ds[
        [v for v in arrangement + others if (v not in drop) and (v in ds.data_vars)]
    ]
    return ds
