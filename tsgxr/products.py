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
        c: f"#{color:06x}"
        for c, color in zip(
            pcls.classes.values(),
            pcls.colors,
        )
    }  #
