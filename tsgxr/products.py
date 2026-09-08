import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray

from .util import Handle

logger = Handle(__name__)


def _product_summary_table(
    ds: xarray.Dataset, which: str, level: str = "Grp"
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
    grps = {
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
    df.name = "sTSA{which}"
    df.columns.name = None
    df.columns.name = next(iter(ds.dims))
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


def bgrint_to_rgb(v: int):
    # fractional RGB from bgr integer
    # https://github.com/AuScope/nvcl_kit/blob/2ab72a9c2133715a1ffc75af282e2824b2681bca/nvcl_kit/reader.py#L62-L68
    return ((v & 255) / 255.0, ((v & 65280) >> 8) / 255.0, (v >> 16) / 255.0)


def get_product_colormap(ds: xarray.Dataset, which: str = "sTSAS"):
    pcls = next(
        iter(
            [
                v
                for ix, v in ds.attrs["class"].items()
                # these seem to be _ delimited
                if f"{which[0]}_{which[1:]}".upper() in v.name.upper()
            ]
        )
    )
    return {
        c: bgrint_to_rgb(color)
        for c, color in zip(
            pcls.classes.values(),
            pcls.colors,
        )
    }


def plot_product_downhole(
    ds: xarray.Dataset,
    which: str,
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
    products = products_to_mineral_table(ds, which=which)
    colormap = get_product_colormap(ds, which=which)
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
