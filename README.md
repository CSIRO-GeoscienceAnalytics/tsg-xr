# `tsg-xr`: A tool for loading TSG datasets into Xarray

The file format associated with [The Spectral Geologist™](https://research.csiro.au/thespectralgeologist/) 
(and specifically [Hylogger™](https://corescan.com.au/products/hylogger/) datasets which
have been processed with the software) consists of an ensemble of files:
* Binary data files containing spectra, high resolutoin imagery and profilometer data
* Configuration files (principally text, similar in format to TOML)
* Low resolution core imagery exports (hole overview, per-tray imagery; as JPEG images with associated markup)

`tsg-xr` heavily leverages the filereader of [`pytsg`](https://https://github.com/Geological-Survey-of-Western-Australia/pytsg) to 
provide access to these data, and presents data in an [Xarray](https://xarray.pydata.org) format to condense the 
otherwise complex arrangement. Here `pytsg` provides an efficient interface to the 
binary components of the TSG file format, and `tsg-xr` is largely just arranging this into a condensed 
data structure which allows easier subseqent use (and serialization to indexable formats, e.g. 
[Zarr](https://zarr.dev/)).

## Usage

`tsg-xr` is intended to be used to read directories containing ensembles of TSG files; 
to do so just point the `open_tsg` funnction at the appropriate directory:

```python
from tsgxr import open_tsg

DT : xarray.DataTree = open_tsg("./Hylogger_Hole_42")
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
spectra : xarray.DataArray = DT["NIR"].ds["Spectra"]

spectra.sel(depth=slice(9.2, 9.3)).plot.line(
    x="wavelength", add_legend=False, color="k", alpha=0.51
)
```

---

Scalars and other spectral features are also available; spectral feature (centre, depth, width) data is grouped for brevity:

```python
DT["NIR"].ds.Centres
DT["NIR"].ds.Depths
DT["NIR"].ds.Widths
DT["NIR"].ds["Grp1 sTSAS"]
...
DT["NIR"].ds["Min1 sTSAS"]
...
DT["NIR"].ds["Wt1 sTSAS"]
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


The package is can be installed from PyPI using `pip` with:
```bash
pip install tsgxr
```

The package is also directly installable from GitHub using `pip` with:
```bash
pip install git+https://github.com/CSIRO-GeoscienceAnalytics/tsg-xr
```

**Option 2: Setup an Environment**

This repository is set up to use `uv` for environment management, and a `uv.lock` file is included in this repository. After cloning this repository and navigating to this directory, `uv venv` followed by `./.venv\Scripts\activate` and `uv sync` should create a `.venv` virtual environment in which to run `tsg-xr` code. Environments can be activated in terminal environments with `.\.venv\Scripts\activate` (on Windows).

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

## Data Structure Overview 

An example of the data structure used is given below, for the `STAVELY_17` hole available from the NVCL (note that the choice of index coordinate - between depth and sample - will affect 
some of this structure's orientation):

```python
DT : xarray.DataTree = open_tsg(
    "./07e4dcac-5216-44a6-9a6b-0c4c1f7ce7d", index_coord="depth", lazy=True, chunks=512
)
```

```
<xarray.DataTree>
Group: /
├── Group: /NIR
│       Dimensions:                    (depth: 19754, feature: 25, wavelength: 531)
│       Coordinates:
│         * depth                      (depth) float32 79kB 0.004111 0.004112 ... 133.4
│           sample                     (depth) uint64 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│           tray                       (depth) uint16 40kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section                    (depth) uint8 20kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section-part               (depth) uint8 20kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section-position           (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           hole                       (depth) object 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│         * feature                    (feature) <U2 200B '1' '2' '3' ... '23' '24' '25'
│         * wavelength                 (wavelength) float64 4kB 380.0 384.0 ... 2.5e+03
│           band                       (wavelength) int64 4kB dask.array<chunksize=(512,), meta=np.ndarray>
│       Data variables: (12/70)
│           HoleID                     (depth) object 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Date                       (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Depth (m)                  (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Tray                       (depth) <U16 1MB dask.array<chunksize=(512,), meta=np.ndarray>
│           Section                    (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Centres                    (feature, depth) float32 2MB dask.array<chunksize=(25, 512), meta=np.ndarray>
│           ...                         ...
│           TraySamp                   (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           colour mod_sat_intens      (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           core_qual                  (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           prof_min                   (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           sec_end_mask               (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Spectra                    (depth, wavelength) float32 42MB dask.array<chunksize=(512, 512), meta=np.ndarray>
│       Attributes: (12/38)
│           core_qual:                 ['Void', 'Rubble', 'Crack', 'Core']
│           TSA704_S Minerals:         ['Opal', 'Dickite', 'Kaolinite-PX', 'Kaolinite...
│           Tray:                      ['0001', '0002', '0003', '0004', '0005', '0006...
│           TSA704_V Groups:           ['MISC-SILICATE', 'CARBONATE', 'SULPHATE', 'OX...
│           HyLogDiag:                 ['wc ws', 'al wc ws', 'wc ws pz', 'al wc ws pz...
│           RockMarks:                 []
│           ...                        ...
│           tsasettings 1:             {'items_full': '17', 'items_sub': '6', 'trains...
│           sclrsets:                  {}
│           domain 0:                  {'name': 'Default', 'samp0': '0', 'samp1': '23...
│           events:                    {}
│           wavelength specs:          {'start': 380.0, 'end': 2500.0, 'unit': 'nm'}
│           batch:                     {'commands': '1', 'name': 'Kahuna,15', 'descri...
├── Group: /TIR
│       Dimensions:                        (depth: 19754, feature: 25, wavelength: 341)
│       Coordinates:
│         * depth                          (depth) float32 79kB 0.004111 ... 133.4
│           sample                         (depth) uint64 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│           tray                           (depth) uint16 40kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section                        (depth) uint8 20kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section-part                   (depth) uint8 20kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section-position               (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           hole                           (depth) object 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│         * feature                        (feature) <U2 200B '1' '2' '3' ... '24' '25'
│         * wavelength                     (wavelength) float64 3kB 6e+03 ... 1.45e+04
│           band                           (wavelength) int64 3kB dask.array<chunksize=(341,), meta=np.ndarray>
│       Data variables: (12/63)
│           HoleID                         (depth) object 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Date                           (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Depth (m)                      (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Tray                           (depth) <U16 1MB dask.array<chunksize=(512,), meta=np.ndarray>
│           Section                        (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Centres                        (feature, depth) float32 2MB dask.array<chunksize=(25, 512), meta=np.ndarray>
│           ...                             ...
│           SecSamp                        (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Subpix                         (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           TIRDeltaTemp                   (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           TirBkgOffset                   (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           TraySamp                       (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           Spectra                        (depth, wavelength) float32 27MB dask.array<chunksize=(512, 341), meta=np.ndarray>
│       Attributes: (12/32)
│           TSA703_T Groups:           ['SILICA', 'K-FELDSPAR', 'PLAGIOCLASE', 'GARNE...
│           TSA703_T Minerals:         ['Opal', 'Quartz', 'Anorthoclase', 'Microcline...
│           HoleID:                    ['stavely_17']
│           Domain:                    ['Default']
│           HyLogDiag:                 ['wc ws', 'al wc ws', 'wc ws pz', 'al wc ws pz...
│           RockMarks:                 []
│           ...                        ...
│           tsatircal:                 {'tirq': '0.000875 0.000890 0.000905 0.000918 ...
│           sclrsets:                  {}
│           domain 0:                  {'name': 'Default', 'samp0': '0', 'samp1': '23...
│           events:                    {}
│           wavelength specs:          {'start': 6000.0, 'end': 14500.0, 'unit': 'nm'}
│           batch:                     {'commands': '12', 'name': 'Restrahlen_feature...
├── Group: /Lidar
│       Dimensions:           (depth: 19754)
│       Coordinates:
│         * depth             (depth) float32 79kB 0.004111 0.004112 ... 133.4 133.4
│           sample            (depth) uint64 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│           tray              (depth) uint16 40kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section           (depth) uint8 20kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section-part      (depth) uint8 20kB dask.array<chunksize=(512,), meta=np.ndarray>
│           section-position  (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
│           hole              (depth) object 158kB dask.array<chunksize=(512,), meta=np.ndarray>
│       Data variables:
│           Lidar             (depth) float32 79kB dask.array<chunksize=(512,), meta=np.ndarray>
└── Group: /Image
        Dimensions:  (depth: 2898500, channel: 3, width: 926)
        Coordinates:
          * depth    (depth) float32 12MB 0.004111 0.004177 0.004243 ... 156.0 156.0
            section  (depth) uint16 6MB dask.array<chunksize=(512,), meta=np.ndarray>
            tray     (depth) uint16 6MB dask.array<chunksize=(512,), meta=np.ndarray>
          * channel  (channel) int64 24B 0 1 2
          * width    (width) float64 7kB -0.03054 -0.03047 -0.0304 ... 0.03047 0.03054
        Data variables:
            Image    (depth, width, channel) uint8 8GB dask.array<chunksize=(512, 512, 3), meta=np.ndarray>
```

## Performance Overview

Some rough performance numbers are given below, comparing `tsg-xr` and `pytsg`, performed with Python 3.13 on Windows using an i7-13850HX (2.10 GHz) reading from a Gen4 NVME.

Note that `tsg-xr` *is generally slower* to load data (it uses `pytsg` for some of the basic loading steps and data classes), 
but it provides a more formatted/annotated data structure, translation of coordinates, and nodata values.
Further, the main benefits of lazy loading are for true color imagery, and principally for memory usage where you're not 
planning to load the whole dataset (at least at once), including where the imagery is of greater size than availabile RAM. 


### Moderately Sized Dataset: `STAVELY_17`

`STAVELY_17` is a hylogger dataset with ID `07e4dcac-5216-44a6-9a6b-0c4c1f7ce7d` in NVCL shown above, CRAS is 296 MB and it has NIR and TIR spectral data 
totalling 155 MB.

Loading the whole dataset *without an image*:
```python
> %timeit open_tsg(hyloggerdir, image=False, lazy=False) # lazy tsg-xr, loaded
952 ms ± 24.9 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

> %timeit open_tsg(hyloggerdir, image=False) # lazy tsg-xr
847 ms ± 4.15 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

> %timeit read_package(hyloggerdir, read_cras_file=False) # pytsg
420 ms ± 17.3 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Loading the dataset *with an image*:

```python
%timeit open_tsg(hyloggerdir, image=True, lazy=False) # tsg-xr
4.21 s ± 539 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

%timeit open_tsg(hyloggerdir, image=True) # lazy tsg-xr
1.03 s ± 51.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

%timeit read_package(hyloggerdir, read_cras_file=True) # pytsg
9.08 s ± 193 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Reading a TSG spectral dataset:

```python
> %timeit xarray.open_dataset(tsgfile, engine="tsg") # lazy tsg-xr
494 ms ± 20.5 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

> %timeit pytsg.parse_tsg.read_tsg_bip_pair(tsgfile, bipfile, "NIR",) # pytsg
185 ms ± 5.58 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

### Large Dataset: `Barnicarndy 1`

`Barnicarndy 1` is a hylogger dataset with ID `81eabec7-0d6d-4e58-8e8a-d0c579b9839` in NVCL, CRAS is 6.4 GB with NIR and TIR spectral data totalling 2.0 GB. When uncompressed, the contained JPEG imagery is over 220GB, and as such would not fit in memory on standard
machines.

Loading the whole dataset *without an image*:

```python
> %timeit open_tsg(hyloggerdir, image=False) # lazy tsg-xr
9.92 s ± 696 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

> %timeit read_package(hyloggerdir, read_cras_file=False) # pytsg
4.59 s ± 161 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Loading the dataset *with an image* (note: won't fit in memory, so `pytsg` metric not 
given here):

```python
%timeit open_tsg(hyloggerdir, image=True) # lazy tsg-xr
10.2 s ± 71.9 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Reading a TSG spectral dataset:

```python
> %timeit xarray.open_dataset(tsgfile, engine="tsg") # lazy tsg-xr
5.19 s ± 118 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

> %timeit pytsg.parse_tsg.read_tsg_bip_pair(tsgfile, bipfile, "NIR",) # pytsg
2.46 s ± 87.3 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```
