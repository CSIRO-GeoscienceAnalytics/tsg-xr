import xarray

from .util import Handle

logger = Handle(__name__)


def _product_summary_table(ds, which, level="Grp"):
    grps = {
        ix + 1: v
        for ix, v in enumerate(
            [
                v.split(" ")[0]
                for v in ds.data_vars
                if (v.startswith(level) and v.endswith(f"sTSA{which}"))
            ]
        )
    }
    df = (
        sum(
            [
                ds[[f"{g} sTSA{which}", f"Wt{ix} sTSA{which}"]]
                .to_dataframe()
                .reset_index(drop=True)
                .pivot(columns=f"{g} sTSA{which}", values=f"Wt{ix} sTSA{which}")
                .fillna(0)
                for ix, g in grps.items()
            ]
        )
    ).set_index(ds.depth.values)
    df.name = "sTSA{which}"
    df.columns.name = None
    return df.where(df > 0).dropna(how="all", axis=1)


def products_to_group_table(ds: xarray.Dataset, which="S"):
    return _product_summary_table(ds, which, level="Grp")


def products_to_mineral_table(ds: xarray.Dataset, which="S"):
    return _product_summary_table(ds, which, level="Min")
