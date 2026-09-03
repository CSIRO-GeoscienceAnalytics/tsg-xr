import struct
from pathlib import Path

import joblib
import numpy as np
import xarray
from pytsg.parse_tsg import CrasHeader, SectionInfo, TrayInfo
from simplejpeg import decode_jpeg

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


class CRASBackend(xarray.backends.BackendEntrypoint):
    """
    A parallel-loading backend to open TSG truecolor imagery.
    """

    description = "Load TSG truecolor imagery using xarray"

    def open_dataset(
        self,
        filename_or_obj,
        header_format="20s2I8h4I2h",
        tray_info_format: str = "3f2i",
        section_info_format: str = "4f3i",
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

            diff_offset = np.diff(self.offsets, prepend=1).astype(np.uint64)
            overflow_finder = np.where(diff_offset < -1)[0].astype(np.uint64)
            if len(overflow_finder) > 1:
                raise IndexError("Chunk offset array wraps around more than once")

            if len(overflow_finder) > 0:
                # add np.int32 max to the offset array this should be ok, unless there is a case where there is more than 1 overflow,
                # in which case I expect the cras reading component to crash
                self.offsets[overflow_finder[0] :] += np.int64(
                    np.iinfo(np.uint32).max + 1
                )

            curpos: int = 0
            # reading the file only takes a small amount of the time
            datas = []
            for i in np.arange(self.header.nchunks):
                img_nbytes = self.offsets[i + 1] - self.offsets[i]
                file.seek(
                    (self.offsets[i] + 4 * (self.header.nchunks + 1) + 64).astype(int)
                )
                datas += [
                    (curpos, file.read(img_nbytes))
                ]  # the pos is only needed for later allocation workflows... not if stacking
                curpos += self.header.chunksize

            # not sure if it's faster to alloate the array upfront or just concatenate the arrays

            # def assign_data(pos, data):
            #     img_data = decode_jpeg(data, colorspace="BGR")[::-1]
            #     cras[pos : (pos + img_data.shape[0]), :, :] = img_data
            #     del img_data

            def get_img(pos, data):
                return decode_jpeg(data)[::-1]

            cras = np.vstack(
                list(
                    joblib.Parallel(
                        n_jobs=-1, backend="threading", return_as="generator"
                    )(
                        (joblib.delayed(get_img)(*d) for d in datas),
                    )
                )
            )
            del datas
            # cras = np.zeros(
            #     (self.header.nl, self.header.ns, self.header.nb), dtype=np.uint8
            # )
            # for i in np.arange(self.header.nchunks):
            #     img_nbytes = self.offsets[i + 1] - self.offsets[i]
            #     file.seek(
            #         (self.offsets[i] + 4 * (self.header.nchunks + 1) + 64).astype(int)
            #     )
            #     img_data = decode_jpeg(file.read(img_nbytes), colorspace="BGR")[::-1]
            #     # reverse the channels, and flip the image upsidedown, after decoding
            #     cras[curpos : (curpos + img_data.shape[0]), :, :] = img_data
            #     curpos += img_data.shape[0]
            #     del img_data
            #     if pbar:
            #         pbar.update(1)

            info_table_start = (
                64
                + (self.header.nchunks + 1) * 4
                + self.offsets[self.header.nchunks]
                - self.offsets[0]
            ).astype(np.uint64)
            file.seek(info_table_start)

            self.tray: list[TrayInfo] = []
            for i in range(self.header.ntrays):
                bytes = file.read(20)
                self.tray.append(TrayInfo(*struct.unpack(tray_info_format, bytes)))

            self.section: list[SectionInfo] = []
            for i in range(self.header.nsections):
                bytes = file.read(28)
                self.section.append(
                    SectionInfo(*struct.unpack(section_info_format, bytes))
                )

        da = xarray.DataArray(
            data=cras,
            dims=("x", "y", "channel"),
            coords={
                "section": (
                    "x",
                    np.hstack(
                        [
                            np.ones(s.nlines, dtype="int32") * ix
                            for ix, s in enumerate(self.section)
                        ]
                    ),
                ),
                "tray": (
                    "x",
                    np.hstack(
                        [
                            np.ones(s.nlines, dtype="int32") * ix
                            for ix, s in enumerate(self.tray)
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
        tray_info_format: str = "3f2i",
        section_info_format: str = "4f3i",
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

            diff_offset = np.diff(self.offsets, prepend=1).astype(np.uint64)
            overflow_finder = np.where(diff_offset < -1)[0].astype(np.uint64)
            if len(overflow_finder) > 1:
                raise IndexError("Chunk offset array wraps around more than once")

            if len(overflow_finder) > 0:
                # add np.int32 max to the offset array this should be ok, unless there is a case where there is more than 1 overflow,
                # in which case I expect the cras reading component to crash
                self.offsets[overflow_finder[0] :] += np.int64(
                    np.iinfo(np.uint32).max + 1
                )
            info_table_start = (
                64
                + (self.header.nchunks + 1) * 4
                + self.offsets[self.header.nchunks]
                - self.offsets[0]
            ).astype(np.uint64)
            file.seek(info_table_start)

            self.tray: list[TrayInfo] = []
            for i in range(self.header.ntrays):
                bytes = file.read(20)
                self.tray.append(TrayInfo(*struct.unpack(tray_info_format, bytes)))

            self.section: list[SectionInfo] = []
            for i in range(self.header.nsections):
                bytes = file.read(28)
                self.section.append(
                    SectionInfo(*struct.unpack(section_info_format, bytes))
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
                "x": np.arange(arr.shape[0]) + self.header.chunksize * chunkidx_start,
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
                            np.ones(s.nlines, dtype="int32") * ix
                            for ix, s in enumerate(backend_array.section)
                        ]
                    ),
                ),
                "tray": (
                    "x",
                    np.hstack(
                        [
                            np.ones(s.nlines, dtype="int32") * ix
                            for ix, s in enumerate(backend_array.tray)
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
