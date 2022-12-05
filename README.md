# `tsg-xr`: A tool for loading TSG datasets into Xarray

The file format associated with The Spectral Geologist™ (and particularly with Hylogger™ datasets which 
have been processed with the software) consists of an ensemble of files:
* Binary data files containing spectra, high resolutoin imagery and profilometer data
* Configuration files (principally text, similar in format to TOML)
* Low resolution core imagery exports (hole overview, per-tray imagery; as JPEG images with associated markup)

`tsg-xr` heavily leverages the filereader of [`pytsg`](https://github.com/FractalGeoAnalytics/pytsg) to 
provide access to these data, and presents data in an [Xarray](xarray.pydata.org) format to condense the 
otherwise complex arrangement and separation of components. `pytsg` provides an efficient interface to the 
binary components of the TSG file format, and `tsg-xr` is largely just arranging this into a condensed 
data structure which allows easier subseqent use (and serialization to indexable formats, e.g. 
[Zarr](https://zarr.readthedocs.io).


## Installation 

The `tsg-xr` pacakge can be installed standalone into your local environment using `pip`, or you can create an 
environment with related dependencies using Anaconda (useful for a development scenario, or if you're only using
the tool for a singular project).

**Option 1: Standalone Installation**

The package is also directly installable from GitHub using `pip` with:
```bash
pip install git+https://github.com/CSIRO-GeoscienceAnalytics/tsg-xr
```

**Option 1: Setup an Environment**

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

