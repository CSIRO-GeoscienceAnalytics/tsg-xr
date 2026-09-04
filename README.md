# `tsg-xr`: A tool for loading TSG datasets into Xarray

The file format associated with [The Spectral Geologist™](https://research.csiro.au/thespectralgeologist/) 
(and specifically [Hylogger™](https://corescan.com.au/products/hylogger/) datasets which
have been processed with the software) consists of an ensemble of files:
* Binary data files containing spectra, high resolutoin imagery and profilometer data
* Configuration files (principally text, similar in format to TOML)
* Low resolution core imagery exports (hole overview, per-tray imagery; as JPEG images with associated markup)

`tsg-xr` heavily leverages the filereader of [`pytsg`](https://https://github.com/Geological-Survey-of-Western-Australia/pytsg) to 
provide access to these data, and presents data in an [Xarray](xarray.pydata.org) format to condense the 
otherwise complex arrangement. Here `pytsg` provides an efficient interface to the 
binary components of the TSG file format, and `tsg-xr` is largely just arranging this into a condensed 
data structure which allows easier subseqent use (and serialization to indexable formats, e.g. 
[Zarr](https://zarr.readthedocs.io)).

## Usage

`tsg-xr` is intended to be used to read directories containing ensembles of TSG files; 
to do so just point the `load_tsg` funnction at the appropriate directory:

```python
from tsgxr import load_tsg

DT : xarray.DataTree = load_tsg("./Hylogger_Hole_42")
```

---

Key array-based data can be accessed directly from this `xarray.Datatree` object:

```python
DT['NIR/Spectra']: xarray.Dataset
DT['NIR/Products']: xarray.Dataset
DT['TIR/Spectra']: xarray.Dataset
DT['TIR/Products']: xarray.Dataset
DT['Image']: xarray.Dataset
DT['Lidar']: xarray.Dataset
```

For example, to extract and plot the first metre of core imagery (note here the 'Image' key is repeated twice due to restrictions on the structure):

```python
import matplotlib.pyplot as plt 

image : xarray.DataArray = DT["Image"].ds.Image
DT["Image"].ds.Image.sel(depth=slice(0, 1)).plot.imshow(yincrease=False)
plt.gca().set(aspect="equal"); # fix the aspect ratio
```

Similarly, to plot the spectra from a specific interval (e.g. 9.2 to 9.3m here) against wavelength, you can provide a slice to the `xarray.DataArray.sel` method:

```python
spectra : xarray.DataArray = DT["NIR/Spectra"].ds["Spectra"]

spectra.sel(depth=slice(9.2, 9.3)).plot.line(
    x="wavelength", add_legend=False, color="k", alpha=0.5
)
```

---

Scalars and other spectral features are also available; spectral feature (centre, depth, width) data is grouped for brevity:

```python
products: xarray.Dataset = DT["NIR"]["Products"].ds

products.Centres
products.Depths
products.Widths
products["Grp1 sTSAS"]
...
products["Min1 sTSAS"]
...
products["Wt1 sTSAS"]
...
```

---

Configuration related to integer-encoding of sample data is also included in the dataset attributes:

```python
products.attrs
```

## Installation 

The `tsg-xr` pacakge can be installed standalone into your local environment using `pip`, or you can create an environment with related dependencies using [`uv`](https://docs.astral.sh/uv/) (useful for a development scenario, or if you're only using the tool for a singular project).

**Option 1: Standalone Installation**

The package is also directly installable from GitHub using `pip` with:
```bash
pip install git+https://github.com/CSIRO-GeoscienceAnalytics/tsg-xr
```

**Option 2: Setup an Environment**

This repository is set up to use `uv` for environment management, and a `uv.lock` file is included in this repository. After cloning this repository and navigating to this directory, `uv venv` followed by `./.venv\Scripts\activate` and `uv sync` should create a `.venv` virtual environment in which to run `tsg-xr` code. Environments can be activated in terminal environments with `.\.venv\Scripts\activate` (on Windows).

## Data Structure Overview 

An example of the data structure used is given below, for the `STAVELY_17` hole available from the NVCL (note that the choice of index coordinate - between depth and sample - will affect 
some of this structure's orientation):

```python
DT : xarray.DataTree = load_tsg(
    "./07e4dcac-5216-44a6-9a6b-0c4c1f7ce7d", index_coord="depth", lazy=True, chunks=512
)
```

```
<xarray.DataTree>
Group: /
├── Group: /NIR
│   ├── Group: /NIR/Spectra
│   │       Dimensions:           (depth: 21350, wavelength: 531)
│   │       Coordinates:
│   │         * depth             (depth) float64 171kB 0.004111 0.004112 ... 156.0 156.0
│   │           sample            (depth) int64 171kB 1 4 2 5 0 ... 23290 23269 23251 23294
│   │           tray              (depth) int64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           section           (depth) int64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           section-part      (depth) int64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           section-position  (depth) float64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           hole              (depth) object 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │         * wavelength        (wavelength) float64 4kB 380.0 384.0 ... 2.496e+03 2.5e+03
│   │       Data variables:
│   │           Spectra           (depth, wavelength) float32 45MB dask.array<chunksize=(508, 512), meta=np.ndarray>
│   └── Group: /NIR/Products
│           Dimensions:           (sample: 23375)
│           Coordinates:
│             * sample            (sample) int64 187kB 0 1 2 3 4 ... 23371 23372 23373 23374
│               tray              (sample) int64 187kB 1 1 1 1 1 1 1 ... 50 50 50 50 50 50
│               section           (sample) int64 187kB 1 1 1 1 1 1 1 1 1 ... 4 4 4 4 4 4 4 4
│               section-part      (sample) int64 187kB 1 2 3 4 5 6 ... 121 122 123 124 125
│               depth             (sample) float64 187kB 0.004117 0.004111 ... 156.0 156.0
│               section-position  (sample) float64 187kB 6.264 14.26 22.26 ... 988.0 996.0
│               hole              (sample) object 187kB 'STAVELY_17' ... 'STAVELY_17'
├── Group: /TIR
│   ├── Group: /TIR/Spectra
│   │       Dimensions:           (depth: 21350, wavelength: 341)
│   │       Coordinates:
│   │         * depth             (depth) float64 171kB 0.004111 0.004112 ... 156.0 156.0
│   │           sample            (depth) int64 171kB 1 4 2 5 0 ... 23290 23269 23251 23294
│   │           tray              (depth) int64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           section           (depth) int64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           section-part      (depth) int64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           section-position  (depth) float64 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │           hole              (depth) object 171kB dask.array<chunksize=(508,), meta=np.ndarray>
│   │         * wavelength        (wavelength) float64 3kB 6e+03 6.025e+03 ... 1.45e+04
│   │       Data variables:
│   │           Spectra           (depth, wavelength) float32 29MB dask.array<chunksize=(508, 341), meta=np.ndarray>
│   └── Group: /TIR/Products
│           Dimensions:           (sample: 23375)
│           Coordinates:
│             * sample            (sample) int64 187kB 0 1 2 3 4 ... 23371 23372 23373 23374
│               tray              (sample) int64 187kB 1 1 1 1 1 1 1 ... 50 50 50 50 50 50
│               section           (sample) int64 187kB 1 1 1 1 1 1 1 1 1 ... 4 4 4 4 4 4 4 4
│               section-part      (sample) int64 187kB 1 2 3 4 5 6 ... 121 122 123 124 125
│               depth             (sample) float64 187kB 0.004117 0.004111 ... 156.0 156.0
│               section-position  (sample) float64 187kB 6.264 14.26 22.26 ... 988.0 996.0
│               hole              (sample) object 187kB 'STAVELY_17' ... 'STAVELY_17'
├── Group: /Lidar
│       Dimensions:           (depth: 21350, sample: 23375)
│       Coordinates:
│           depth             (sample) float64 187kB 0.004117 0.004111 ... 156.0 156.0
│         * sample            (sample) int64 187kB 0 1 2 3 4 ... 23371 23372 23373 23374
│           tray              (sample) int64 187kB 1 1 1 1 1 1 1 ... 50 50 50 50 50 50
│           section           (sample) int64 187kB 1 1 1 1 1 1 1 1 1 ... 4 4 4 4 4 4 4 4
│           section-part      (sample) int64 187kB 1 2 3 4 5 6 ... 121 122 123 124 125
│           section-position  (sample) float64 187kB 6.264 14.26 22.26 ... 988.0 996.0
│           hole              (sample) object 187kB 'STAVELY_17' ... 'STAVELY_17'
│       Data variables:
│           Lidar             (depth) float32 85kB 92.44 75.37 75.98 ... 58.43 1.3
└── Group: /Image
        Dimensions:  (depth: 2898500, width: 926, channel: 3, x: 2898500)
        Coordinates:
          * depth    (depth) float32 12MB 0.004111 0.004177 0.004243 ... 156.0 156.0
          * width    (width) float64 7kB -0.03054 -0.03047 -0.0304 ... 0.03047 0.03054
          * channel  (channel) int64 24B 0 1 2
            section  (x) int32 12MB dask.array<chunksize=(512,), meta=np.ndarray>
            tray     (x) int32 12MB dask.array<chunksize=(512,), meta=np.ndarray>
        Dimensions without coordinates: x
        Data variables:
            Image    (depth, width, channel) uint8 8GB dask.array<chunksize=(511, 512, 3), meta=np.ndarray>
```

## Command Line Interface

### Converting TSG files to Zarr

A minimal command line interface exists for converstion of TSG files to Zarr archives. 
Generally, if you're using `uv`, you would use `uv run` to do this; otherwise you can directly use the `tsgxr` entry point where the respective environment is activated.

A selection of configuration options are avialable from the commandline, which can be found under the help menu:

```bash
uv run tsgxr tsg2zarr --help
```
Basic usage is as follows, where `<Path>` refers to either i) an individual TSG scalars file (`.tsg`), ii) a Hylogger TSG directory, or iii) a directory containing multiple Hylogger TSG directories (multiple datasets can be converted simultaneously):

```bash
uv run tsgxr tsg2zarr <Path>
```

Outputs are by default added to the Hylogger TSG directories themselves, but can be optionally collated into a separate directory; outputs will use the hole name extracted from the TSG dataset and be specific to the spectra specified (NIR or TIR):

```bash
uv run tsgxr tsg2zarr <Path> --output_dir "./collated_zarr_archives/"
```

Note that by default, this will create zipped Zarr archives. These can be directly opened in e.g. Xarray.