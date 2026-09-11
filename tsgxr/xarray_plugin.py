import struct
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import xarray
from pytsg.parse_tsg import (
    CrasHeader,
    _calculate_wavelengths,
    _find_header_sections,
    _parse_tsg,
    _read_tsg_file,
)
from simplejpeg import decode_jpeg

from .read import coords_from_sampleheaders, product_dataset_to_xarray
from .util import Handle

logger = Handle(__name__)


try:
    import dask

    def get_lock():
        return dask.utils.SerializableLock()

except ImportError:
    # TODO: Could use threading.lock on linux, there are alternates on windows
    from contextlib import contextmanager

    def get_lock():
        @contextmanager
        def mgr():
            try:
                yield None
            finally:
                pass

        return mgr()


def parse_scalars(
    scalars: np.ndarray,
    classes: "list[ClassHeaders]",
    bandheaders: "list[BandHeaders]",
    nodata: int = -1,
) -> pd.DataFrame:
    """
    Map scalar values to classes, where appropriate.
    """
    df = pd.DataFrame(
        {
            band.name: (
                np.where(
                    np.isclose(scalars[:, band.band], np.finfo("float32").min),
                    -1,
                    scalars[:, band.band],
                )
            )
            for band in bandheaders
        }
    )

    for band in bandheaders:
        if (band.flag == 2) & (band.class_number > -1):
            df[band.name] = pd.Series(df[band.name].astype(np.int16)).map(
                classes[int(band.class_number)].classes
            )
    return df


class BIPBackendArray(xarray.backends.BackendArray):
    def __init__(
        self,
        filename_or_obj,
        shape=None,
        dtype=None,
        lock=None,
        chunks=None,
        n_bands: int = 512,
        n_samples: int | None = None,
    ):

        self.lock = lock
        if chunks is None:
            chunks = {}

        fpath = Path(filename_or_obj)
        self.tsg, self.bip = None, None
        self.bip = fpath if fpath.suffix == ".bip" else fpath.with_suffix(".bip")
        self.tsg = fpath if fpath.suffix == ".tsg" else fpath.with_suffix(".tsg")
        if not self.bip.exists() and self.tsg.exists():
            raise FileNotFoundError(
                f"Missing file: {','.join(([self.bip.name] if not self.bip.exists() else []) + ([self.tsg.name] if not self.tsg.exists() else []))}"
            )
        self.fstr = _read_tsg_file(self.tsg)
        self.headers = _find_header_sections(self.fstr)
        self.info = _parse_tsg(self.fstr, self.headers)
        self.wavelength = _calculate_wavelengths(
            self.info["wavelength specs"], self.info["coordinates"]
        )
        self.info["coordinates"] = {
            k: int(v) for k, v in self.info["coordinates"].items()
        }
        self.coords = coords_from_sampleheaders(
            self.info["sample headers"], self.wavelength
        )
        self.shape = (
            2,
            self.info["coordinates"]["lastsample"],
            self.info["coordinates"]["lastband"],
        )
        self.dtype = np.dtype(np.float32)

    def __getitem__(self, key: tuple):
        return xarray.core.indexing.explicit_indexing_adapter(
            key,
            self.shape,
            xarray.core.indexing.IndexingSupport.BASIC,
            self._raw_indexing_method,
        )

    def _raw_indexing_method(self, key: tuple):
        key1 = key[1]
        if isinstance(key1, slice):
            start = key1.start or 0
            stop = key1.stop or self.shape[1]  # final pixel
        else:
            start = key1
            stop = key1

        nsamples = stop - start

        with self.lock, open(self.bip, "rb") as f:
            arr = np.fromfile(
                f,
                offset=start
                * self.info["coordinates"]["lastband"]
                * 2
                * np.dtype(self.dtype).itemsize,
                count=nsamples * self.info["coordinates"]["lastband"] * 2,
                dtype=self.dtype,
            ).reshape(2, nsamples, self.info["coordinates"]["lastband"])

        # these need to be integer-indexed
        arr = xarray.DataArray(
            arr,
            dims=("half", "sample", "wavelength"),
            coords={
                "half": np.arange(2),
                "sample": np.arange(start, stop, dtype="uint64"),
                "wavelength": np.arange(self.wavelength.size),
            },
        )
        if isinstance(key, int):
            arr = arr.squeeze()
        return arr.loc[*key].values


class TSGBIPBackend(xarray.backends.BackendEntrypoint):
    """
    A lazy-loading xarray backend to open a single TSG spectral dataset.
    """

    description = "Lazy-load TSG spectral datasets using xarray"

    def open_dataset(
        self,
        filename_or_obj,
        header_format="20s2I8h4I2h",
        tray_info_format: str = "3f2i",
        section_info_format: str = "4f3i",
        drop_variables=None,
        lock=None,
        chunks=None,
    ) -> xarray.Dataset:
        if chunks is None:
            chunks = {}
        backend_array = BIPBackendArray(
            filename_or_obj=filename_or_obj, lock=get_lock()
        )
        self.lock = lock or get_lock()

        # lazy data array representing the spectral array
        da = xarray.DataArray(
            data=xarray.core.indexing.LazilyIndexedArray(backend_array),
            dims=("half", "sample", "wavelength"),
            coords={
                **backend_array.coords,
                "half": np.arange(2),
                "sample": np.arange(
                    0, backend_array.info["coordinates"]["lastsample"], dtype="uint64"
                ),
                "wavelength": backend_array.wavelength,
                "band": ("wavelength", np.arange(backend_array.wavelength.size)),
            },
        )

        product_data = product_dataset_to_xarray(  # products always loads
            parse_scalars(
                da[1].values,
                backend_array.info["class"],
                backend_array.info["band headers"],
            ),
            backend_array.info["class"],
        )
        # handle nodata in spectra here # TODO: can we pull this from the mask?
        ds = product_data.assign(
            Spectra=da[0].where(lambda x: ~np.isclose(x, np.finfo(np.float32).min))
        )
        # all of the useful information from the band headers, sample headers,
        # class headers is incorporated already
        ds.attrs.update(
            {
                k: v
                for k, v in backend_array.info.items()
                if k not in ["band headers", "class", "sample headers"]
            }
        )
        return ds

    def guess_can_open(self, filename_or_obj: str | Path) -> bool:
        fpath = Path(filename_or_obj)
        return ((fpath.suffix == ".bip") and ("tsg" in fpath.stem)) or (
            (fpath.suffix == ".tsg") and ("tsg" in fpath.stem)
        )


class CRASBackend(xarray.backends.BackendEntrypoint):
    """
    A parallel-loading backend to open TSG truecolor imagery.
    """

    description = "Load TSG truecolor imagery using xarray"

    def open_dataset(
        self,
        filename_or_obj,
        header_format="20s2I8h4I2h",
        drop_variables=None,
        lock=None,
    ) -> xarray.Dataset:
        self.filename_or_obj = filename_or_obj
        self.lock = lock or get_lock()
        with self.lock, open(self.filename_or_obj, "rb") as file:
            self.header = CrasHeader(*struct.unpack(header_format, file.read(64)))
            file.seek(64)
            self.offsets = np.ndarray(
                (self.header.nchunks + 1),
                np.uint32,
                file.read(4 * (self.header.nchunks + 1)),
            ).astype(np.uint64)  # deal with +4gb cras files by using uint64

            diff_offset = np.diff(self.offsets, prepend=1)  # this needs to be signed
            overflow_finder = np.where(diff_offset < -1)[0].astype(np.uint64)
            if overflow_finder.size > 1:
                raise IndexError("Chunk offset array wraps around more than once")

            if overflow_finder.size:
                # add np.int32 max to the offset array this should be ok, unless there is a case where there is more than 1 overflow,
                # in which case I expect the cras reading component to crash
                self.offsets[overflow_finder[0] :] += np.uint64(
                    np.iinfo(np.uint32).max + 1
                )

            # reading the file only takes a small amount of the time
            datas = []
            for i in np.arange(self.header.nchunks):
                file.seek(
                    (self.offsets[i] + 4 * (self.header.nchunks + 1) + 64).astype(int)
                )
                datas.append(file.read(self.offsets[i + 1] - self.offsets[i]))

            # not sure if it's faster to alloate the array upfront or just concatenate the arrays
            def get_img(data):
                return decode_jpeg(data, colorspace="BGR")[::-1]

            cras = np.vstack(
                list(
                    joblib.Parallel(
                        n_jobs=-1, backend="threading", return_as="generator"
                    )(
                        (joblib.delayed(get_img)(d) for d in datas),
                    )
                )
            )
            del datas

            info_table_start = (
                64
                + (self.header.nchunks + 1) * 4
                + self.offsets[self.header.nchunks]
                - self.offsets[0]
            ).astype(np.uint64)
            file.seek(info_table_start)

            self.tray = np.rec.array(
                np.fromfile(
                    file,
                    dtype=[
                        ("utlengthmm", "f4"),
                        ("baseheightmm", "f4"),
                        ("coreheightmm", "f4"),
                        ("nsections", "u4"),
                        ("nlines", "u4"),
                    ],
                    count=self.header.ntrays,
                )
            )

            self.section = np.rec.array(
                np.fromfile(
                    file,
                    dtype=[
                        ("utlengthmm", "f4"),
                        ("startmm", "f4"),
                        ("endmm", "f4"),
                        ("trimwidthmm", "f4"),
                        ("startcol", "u4"),
                        ("endcol", "u4"),
                        ("nlines", "u4"),
                    ],
                    count=self.header.nsections,
                )
            )

        da = xarray.DataArray(
            data=cras,
            dims=("x", "y", "channel"),
            coords={
                "section": (
                    "x",
                    np.hstack(
                        [
                            np.ones(n, dtype="uint16") * ix
                            for ix, n in enumerate(self.tray.nlines)
                        ]
                    ),
                ),
                "tray": (
                    "x",
                    np.hstack(
                        [
                            np.ones(n, dtype="uint16") * ix
                            for ix, n in enumerate(self.section.nlines)
                        ]
                    ),
                ),
                "channel": np.arange(3),
            },
        )
        da.attrs.update(
            {"tray": self.tray, "section": self.section},
        )
        return xarray.Dataset({"Image": da})

    def guess_can_open(self, filename_or_obj: str | Path) -> bool:
        return Path(filename_or_obj).name.endswith("cras.bip")


class CRASBackendArray(xarray.backends.BackendArray):
    def __init__(
        self,
        filename_or_obj,
        shape=None,
        dtype=None,
        lock=None,
        chunks=None,
        header_format: str = "20s2I8h4I2h",
    ):
        self.filename_or_obj = filename_or_obj
        self.lock = lock
        if chunks is None:
            chunks = {}
        with self.lock, open(self.filename_or_obj, "rb") as file:
            self.header = CrasHeader(*struct.unpack(header_format, file.read(64)))
            file.seek(64)
            self.offsets = np.ndarray(
                (self.header.nchunks + 1),
                np.uint32,
                file.read(4 * (self.header.nchunks + 1)),
            ).astype(np.uint64)  # deal with +4gb cras files by using uint64

            diff_offset = np.diff(self.offsets, prepend=1)  # this needs to be signed
            overflow_finder = np.where(diff_offset < -1)[0].astype(np.uint64)
            if overflow_finder.size > 1:
                raise IndexError("Chunk offset array wraps around more than once")

            if overflow_finder.size:
                # add np.int32 max to the offset array this should be ok, unless there is a case where there is more than 1 overflow,
                # in which case I expect the cras reading component to crash
                self.offsets[overflow_finder[0] :] += np.uint64(
                    np.iinfo(np.uint32).max + 1
                )

            info_table_start = (
                64
                + (self.header.nchunks + 1) * 4
                + self.offsets[self.header.nchunks]
                - self.offsets[0]
            ).astype(np.uint64)
            file.seek(info_table_start)

            self.tray = np.rec.array(
                np.fromfile(
                    file,
                    dtype=[
                        ("utlengthmm", "f4"),
                        ("baseheightmm", "f4"),
                        ("coreheightmm", "f4"),
                        ("nsections", "u4"),
                        ("nlines", "u4"),
                    ],
                    count=self.header.ntrays,
                )
            )

            self.section = np.rec.array(
                np.fromfile(
                    file,
                    dtype=[
                        ("utlengthmm", "f4"),
                        ("startmm", "f4"),
                        ("endmm", "f4"),
                        ("trimwidthmm", "f4"),
                        ("startcol", "u4"),
                        ("endcol", "u4"),
                        ("nlines", "u4"),
                    ],
                    count=self.header.nsections,
                )
            )

        self.shape = (self.header.nl, self.header.ns, self.header.nb)
        self.imgshape = (self.header.chunksize, self.header.ns, self.header.nb)
        self.dtype = np.dtype(np.uint8)

    def __getitem__(self, key: tuple):
        return xarray.core.indexing.explicit_indexing_adapter(
            key,
            self.shape,
            xarray.core.indexing.IndexingSupport.BASIC,
            self._raw_indexing_method,
        )

    def _raw_indexing_method(self, key: tuple):
        key0 = key[0]
        size = np.dtype(self.dtype).itemsize
        if isinstance(key0, slice):
            start = key0.start or 0
            stop = key0.stop or self.shape[0]  # final pixel
        else:
            start = key0
            stop = key0

        offset = size * start
        # find the relevant chunk
        chunkidx_start = start // self.header.chunksize
        chunkidx_end = stop // self.header.chunksize
        if (stop / self.header.chunksize) > 1.0:
            chunkidx_end += 1
        arrays = []
        with self.lock, open(self.filename_or_obj, "rb") as f:
            for chunk_ID in (
                range(chunkidx_start, chunkidx_end)
                if chunkidx_start != chunkidx_end
                else [chunkidx_start]
            ):
                offset = self.offsets[chunk_ID]
                position = (offset + 4 * (self.header.nchunks + 1) + 64).astype(int)
                length = self.offsets[chunk_ID + 1] - self.offsets[chunk_ID]
                f.seek(position)  # start position
                chunkdata = f.read(length)
                arrays += [np.flipud(decode_jpeg(chunkdata, colorspace="BGR"))]
        arr = np.vstack(arrays)
        # could modify key to offset the start of key[0] by -(self.header.chunksize * chunkidx_start)
        # this is more verbose but useful for debugging where needed...
        arr = xarray.DataArray(
            arr,
            dims=("x", "y", "channel"),
            coords={
                "x": np.arange(arr.shape[0], dtype="uint64")
                + self.header.chunksize * chunkidx_start,
            },
        )
        if isinstance(key, int):
            arr = arr.squeeze()

        return arr.loc[*key].values


class LazyCRASBackend(xarray.backends.BackendEntrypoint):
    """
    A lazy loading backend to open TSG truecolor imagery.
    """

    description = "Open and lazy-load TSG truecolor imagery using xarray"

    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
        dtype=np.int8,
        chunks=None,
    ) -> xarray.Dataset:
        if chunks is None:
            chunks = {}
        backend_array = CRASBackendArray(
            filename_or_obj=filename_or_obj, lock=get_lock()
        )
        da = xarray.DataArray(
            data=xarray.core.indexing.LazilyIndexedArray(backend_array),
            dims=("x", "y", "channel"),
            coords={
                "section": (
                    "x",
                    np.hstack(
                        [
                            np.ones(n, dtype="uint16") * ix
                            for ix, n in enumerate(backend_array.tray.nlines)
                        ]
                    ),
                ),
                "tray": (
                    "x",
                    np.hstack(
                        [
                            np.ones(n, dtype="uint16") * ix
                            for ix, n in enumerate(backend_array.section.nlines)
                        ]
                    ),
                ),
                "channel": np.arange(3),
            },
        )
        da.attrs.update(
            {"tray": backend_array.tray, "section": backend_array.section},
        )
        return xarray.Dataset({"Image": da})

    def guess_can_open(self, filename_or_obj: str | Path) -> bool:
        return Path(filename_or_obj).name.endswith("cras.bip")
