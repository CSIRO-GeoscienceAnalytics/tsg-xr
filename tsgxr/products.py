import re

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray

from .util import Handle

logger = Handle(__name__)


def get_available_systems(ds: xarray.Dataset) -> set:
    """
    Query the set of mineral/group classification systems available on a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to query.

    Returns
    -------
    set
    """
    return {v.split(" ")[1] for v in ds.data_vars if "Grp1" in v or "Min1" in v}


def get_system_subset_attrs(ds: xarray.Dataset, which: str, level=None) -> tuple:
    """
    Get the names of system-related attributes on a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to get attributes from.
    which : str
        Which system to get attributes for. e.g. one of
        'uTSAS', 'uTSAV', 'sTSAS', 'sTSAV', 'uTSAT', 'sTSAT'.
    level : str
        The level to subset attibutes to (either 'Mineral' or 'Groups').

    Returns
    -------
    tuple
    """
    key0 = re.compile(rf"{which[1:-1].upper()}[0-9]+_{which[-1].upper()}")
    # TODO: check if there's independent user and system variables, or just one set?
    items = [
        k
        for k in ds.attrs
        if (
            re.match(key0, k)  # [s/u]TSAS -> TSA704_S
            or (  # e.g. U_SWIR_TSA705
                k.upper().startswith(which[0].upper())
                and (
                    ("SWIR_TSA" in k.upper())
                    if "TSAS" in which
                    else (
                        ("VNIR_TSA" in k.upper()) if "TSAV" in which.upper() else True
                    )
                )
            )
        )
    ]
    if level is not None:
        items = [
            k
            for k in items
            if f"{'Groups' if level.upper().startswith('G') else 'Minerals'}".upper()
            in k.upper()
        ]
    return tuple(items)


def _product_summary_table(
    ds: xarray.Dataset,
    which: str,
    level: str = "Grp",
) -> pd.DataFrame:
    """
    Summarize a TSG scalar/product table, aggregating the long-form
    used in TSG to a full table.

    Parameters
    ----------
    ds : xarray.Dataset
        Product dataset, as loaded in tsg-xr.
    which : str
        Which subset to look at (e.g. S or V for NIR, T for TIR).
    level : str
        Whether to summarize at 'Grp' or 'Min' level.

    Returns
    -------
    pandas.DataFrame
        Dataframe with minerals or groups as columns.
    """
    level = "Grp" if level.upper().startswith("G") else "Min"
    grps = {  # get the groups which correspond to 'which' and 'level'
        ix + 1: v
        for ix, v in enumerate(
            [
                v.split(" ")[0]
                for v in ds.data_vars
                if (v.startswith(level) and v.endswith(f"{which}"))
            ]
        )
    }

    def _get_wideform(ix, g):
        return (
            ds[[f"{g} {which}", f"Wt{ix} {which}"]]
            .to_dataframe()
            .reset_index(drop=True)
            .pivot(columns=f"{g} {which}", values=f"Wt{ix} {which}")
            .fillna(0)
        )

    df = (sum([_get_wideform(ix, g) for ix, g in grps.items()])).set_index(
        ds.depth.values if "depth" in ds.indexes else ds.sample.values
    )
    class_key = next(
        iter(
            [
                k
                for k in get_system_subset_attrs(ds, which=which, level=level)
                if "Colors" not in k
            ]
        )
    )
    df = df[
        [c for c in ds.attrs[class_key] if (c in df.columns)]
    ]  # sort order of columns
    df.name = f"sTSA{which}{'Groups' if level == 'Grp' else 'Minerals'}"
    df.columns.name = None
    df.index.name = next(iter(ds.dims))
    df.columns.name = "{which}group" if level == "Grp" else "{which}mineral"
    return df.where(df > 0).dropna(how="all", axis=1)


def products_to_group_table(ds: xarray.Dataset, which: str = "sTSAS") -> pd.DataFrame:
    """
    Summarize a TSG scalar/product table, aggregating the long-form
    used in TSG to a full table.

    Parameters
    ----------
    ds : xarray.Dataset
        Product dataset, as loaded in tsg-xr.
    which : str
        Which subset to look at (e.g. sTSAS or sTSAV for NIR, T for TIR).

    Returns
    -------
    pandas.DataFrame
        Dataframe with groups as columns.
    """
    return _product_summary_table(ds, which, level="Grp")


def products_to_mineral_table(ds: xarray.Dataset, which: str = "sTSAS") -> pd.DataFrame:
    """
    Summarize a TSG scalar/product table, aggregating the long-form
    used in TSG to a full table.

    Parameters
    ----------
    ds : xarray.Dataset
        Product dataset, as loaded in tsg-xr.
    which : str
        Which subset to look at (e.g. sTSAS or sTSAV for NIR, T for TIR).

    Returns
    -------
    pandas.DataFrame
        Dataframe with minerals as columns.
    """
    return _product_summary_table(ds, which, level="Min")


def get_product_colormap(
    ds: xarray.Dataset,
    which: str = "sTSAS",
    level="Grp",
) -> dict[str, str]:
    """
    Get the product colormap from a dataset for a specific system and level.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to get the colormap from.
    which : str
        Which system to get the colormap for.
    level : str
        Which level to get the colormap for.

    Returns
    -------
    dict [ str, str ]
        Colormap mapping names of minerals/groups to colors (as hex codes,
        as used in `matplotlib`).
    """
    return ds.attrs[
        next(
            iter(
                [
                    k
                    for k in get_system_subset_attrs(ds, which=which, level=level)
                    if "Colors" in k
                ]
            )
        )
    ]


def unstack_arbitrary_feature_data(ds: xarray.Dataset, dtype=np.half) -> xarray.Dataset:
    """
    Unstack the arbitrary spectral feature data from a dataset
    (Centres, Depths, Widths) into a full image dataset with
    dimensions (sample, wavelength).

    Parameters
    ----------
    ds : xarray.Dataset
        Spectral dataset will Centres, Widths and Depths as data variables.
    dtype : type
        Data type to use for the output.

    Returns
    -------
    xarray.Dataset

    Todo
    ----
    * Consider if there's an efficient intermediate form other than the sparse image
      with meaningful coordinates - e.g. long-form features multi-indexed based on depth,
      wavelength.
    """
    idx = "depth" if ("depth" in ds.dims) else "sample"
    complement = "sample" if idx == "depth" else "depth"
    extra_coords = [
        complement,
        "tray",
        "section-part",
        "section-position",
        "hole",
        "section",
    ]
    # note: not all depths or wavelengths might be represented here
    # we need to reindex coordinates along depth/sample as as result
    return (
        ds[["Centres", "Depths", "Widths"]]
        .stack(z=["feature", idx])
        .drop_vars(extra_coords)
        .to_dataframe()
        .reset_index()
        .drop(columns="feature")
        .rename(columns={"Centres": "wavelength"})
        .sort_values([idx, "wavelength"])
        .drop_duplicates(subset=[idx, "wavelength"])
        .set_index([idx, "wavelength"])
        .astype(dtype)
        .to_xarray()
    ).assign_coords(
        {
            c: (idx, ds.sel({idx: ds[idx]}, method="nearest")[c].values)
            if idx in ds[c].dims
            else ds[c]
            for c in extra_coords
            if c in ds
        }
    )


def plot_product_downhole(
    ds: xarray.Dataset,
    which: str,
    level: str = "Grp",
    step: float | None = None,
    ax: matplotlib.axes.Axes | None = None,
    invert: bool = True,
):
    """
    Make a downhole plot of spectral product/scalar data.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing products. Whichever coordinate is used
        as an index here will be used in the downhole dimension.
    which : str
        Which product set to plot (e.g. 'sjCLST').
    level : str
        Whether to summarize at 'Grp' or 'Min' level.
    step : float
        Step to use for compositing, if any. In metres where the
        dataset supplied is indexed by depth (any step > 0.02 makes sense),
        else in sample numbers (i.e. you want to use a step >= 2).
    ax : matplotlib.axes.Axes
        Existing axis to plot on, if one already exists.
    invert : bool
        Whether to invert the yaxis on the plot such that depth
        increases downwards.

    Returns
    -------
    matplotlib.axes.Axes
    """
    products = _product_summary_table(ds, which=which, level=level)
    colormap = get_product_colormap(ds, which=which, level=level)
    if step is not None:
        nsteps = (products.index.values[-1] - products.index.values[0]) // step + 1
        products = products.groupby(
            pd.cut(products.index, products.index.values[0] + np.arange(nsteps) * step),
            observed=False,
        ).agg(lambda x: np.nansum(x) / len(x))
        products.index = products.index.map(
            dict(
                zip(
                    products.index.categories,
                    products.index.categories.map(lambda c: (c.left + c.right) / 2),
                )
            )
        ).values

    if ax is None:
        _fig, ax = plt.subplots(1, figsize=(6, 15))

    bottom = 0
    for c in products.columns:
        ax.barh(
            y=products.index,
            width=products[c],
            label=c,
            height=step,
            left=bottom,
            color=colormap[c],
        )
        bottom += products[c].values

    ax.legend(bbox_to_anchor=(1, 1), loc="upper left", frameon=False)
    ax.set(ylabel="Depth" if next(iter(ds.dims)) == "depth" else "Sample")
    ax.set(title=which)
    if invert:
        ax.invert_yaxis()
    return ax
