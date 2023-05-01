# `tsg-xr`: A tool for loading TSG datasets into Xarray

The file format associated with [The Spectral Geologist™](https://research.csiro.au/thespectralgeologist/) 
(and specifically [Hylogger™](https://corescan.com.au/products/hylogger/) datasets which
have been processed with the software) consists of an ensemble of files:
* Binary data files containing spectra, high resolutoin imagery and profilometer data
* Configuration files (principally text, similar in format to TOML)
* Low resolution core imagery exports (hole overview, per-tray imagery; as JPEG images with associated markup)

`tsg-xr` heavily leverages the filereader of [`pytsg`](https://github.com/FractalGeoAnalytics/pytsg) to 
provide access to these data, and presents data in an [Xarray](xarray.pydata.org) format to condense the 
otherwise complex arrangement. Here `pytsg` provides an efficient interface to the 
binary components of the TSG file format, and `tsg-xr` is largely just arranging this into a condensed 
data structure which allows easier subseqent use (and serialization to indexable formats, e.g. 
[Zarr](https://zarr.readthedocs.io)).

## Usage

`tsg-xr` is intended to be use to read directories containing ensembles of TSG files; to do so just point the `load_tsg` funnction at the appropriate directory:
```python
from tsgxr import load_tsg

ds = load_tsg("./Hylogger_Hole_42")
```
---
Key array-based data can be accessed directly from this `xarray.Dataset` object:
```python
ds.Spectra
ds.Image
ds.Lidar
```

For example, to extract and plot the first metre of core imagery:
```python
import matplotlib.pyplot as plt 								
ds.Image.sel(depth=slice(0, 1)).plot.imshow(yincrease=False)
plt.gca().set(aspect="equal"); # fix the aspect ratio
```

Similarly, to plot the spectra from a specific interval (e.g. 9.2 to 9.3m here) against wavelength, you can provide a slice to the `xarray.DataArray.sel` method (note here the `holedepth` coordinate which is associated with spectral samples, as opposed to the `depth` coordinate assocaited with RGB imagery - they are thus far separate indexes):
```python
ds.Spectra.sel(holedepth=slice(9.2, 9.3)).plot.scatter(x='wavelength', add_legend=False, color='k', alpha=0.5, s=2)
```

---

Scalars and other spectral features are also available; spectral feature (centre, depth, width) data is grouped 
for brevity:
```python
ds.Centres
ds.Depths
ds.Widths
ds["Grp1 sTSAS"]
...
ds["Min1 sTSAS"]
...
ds["Wt1 sTSAS"]
...
```
---
Configuration related to integer-encoding of sample data is also included in the dataset attributes:
```python
ds.attrs
````

## Installation 

The `tsg-xr` pacakge can be installed standalone into your local environment using `pip`, or you can create an 
environment with related dependencies using Anaconda (useful for a development scenario, or if you're only using
the tool for a singular project).

**Option 1: Standalone Installation**

The package is also directly installable from GitHub using `pip` with:
```bash
pip install git+https://github.com/CSIRO-GeoscienceAnalytics/tsg-xr
```

**Option 2: Setup an Environment**

An `environment.yml` file is included in this repository, allowing the creation of a `conda` environment 
where an Anaconda distribution of some form is used. After cloning this repository and navigating to this 
directory, the environment can be created as follows:

```bash
conda env create -f environment.yml
```

Alternatively, if you have `mamba` installed locally (encouraged), you can get there faster with:
```bash
mamba env create -f environment.yml
```

## Command Line Interface

### Converting TSG files to Zarr

A minimal command line interface exists for converstion of TSG files to Zarr archives. A selection of configuration options are avialable from the commandline, which can be found under the help menu:
```bash
tsgxr tsg2zarr --help
```
Basic usage is as follows, where `<Path>` refers to either i) an individual TSG scalars file (`.tsg`), ii) a Hylogger TSG directory, or iii) a directory containing multiple Hylogger TSG directories (multiple datasets can be converted simultaneously):
```bash
tsgxr tsg2zarr <Path>
```
Outputs are by default added to the Hylogger TSG directories themselves, but can be optionally collated into a separate directory; outputs will use the hole name extracted from the TSG dataset and be specific to the spectra specified (NIR or TIR):
```bash
tsgxr tsg2zarr <Path> --output_dir "./collated_zarr_archives/"
```
Note that by default, this will create zipped Zarr archives. These can be directly opened in e.g. Xarray.