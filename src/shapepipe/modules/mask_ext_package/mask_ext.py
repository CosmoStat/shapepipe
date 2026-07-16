"""MASK EXT.

This module contains a class to rasterize an external healsparse mask onto an
image pixel grid, producing the ShapePipe per-image pixel flag file.

:Author: Cail Daley

"""

import re

import numpy as np
from astropy import wcs

from shapepipe.pipeline import file_io

# Per-chunk pixel budget: chunk the image in row bands so that no more than
# this many pixel coordinates are held / queried at once. A tile is ~1e8 px and
# an exposure CCD ~1e7 px; this bound keeps peak memory to a single band.
DEFAULT_PIXEL_BUDGET = 10_000_000


class MaskExt(object):
    """Mask Ext.

    Rasterize an external healsparse mask onto an image pixel grid and write
    the resulting ShapePipe flag file.

    Parameters
    ----------
    image_path : str
        Path to the image whose pixel grid and WCS define the output
    mask_path : str
        Path to the external healsparse mask file
    bit_flag_map : dict
        Mapping from healsparse pixel bit value (int) to output flag value
        (int); a pixel carrying several bits gets the bitwise-OR of the mapped
        flags
    image_num : str
        File number string, inserted into the output file name
    output_dir : str
        Path to the output directory
    w_log : logging.Logger
        Log file
    off_map_flag : int, optional
        Flag value for pixels outside the healsparse footprint (sentinel
        pixels); default is ``0``
    path_external_flag : str, optional
        Path to an external instrument flag file to sum into the mask; default
        is ``None`` (not used)
    image_prefix : str, optional
        Prefix prepended to the output file name base ``flag``; specify
        ``'none'`` or ``''`` for no prefix, default is ``''``
    outname_base : str, optional
        Output file name base, default is ``flag``
    chunk_size : int, optional
        Number of image rows evaluated per chunk; default is derived from the
        image width and :data:`DEFAULT_PIXEL_BUDGET`
    hdu : int, optional
        HDU of the external instrument flag FITS file; default is ``0``

    """

    def __init__(
        self,
        image_path,
        mask_path,
        bit_flag_map,
        image_num,
        output_dir,
        w_log,
        off_map_flag=0,
        path_external_flag=None,
        image_prefix="",
        outname_base="flag",
        chunk_size=None,
        hdu=0,
    ):
        self._image_path = image_path
        self._mask_path = mask_path
        self._bit_flag_map = bit_flag_map
        self._img_number = image_num
        self._output_dir = output_dir
        self._w_log = w_log
        self._off_map_flag = int(off_map_flag)
        self._path_external_flag = path_external_flag

        if (image_prefix.lower() != "none") and (image_prefix != ""):
            self._img_prefix = f"{image_prefix}_"
        else:
            self._img_prefix = ""

        self._outname_base = outname_base
        self._chunk_size = chunk_size
        self._hdu = hdu

        self._set_image_coordinates()

    @staticmethod
    def parse_bit_flag_map(map_string):
        """Parse Bit Flag Map.

        Parse the ``BIT_FLAG_MAP`` config string into a dictionary.

        Parameters
        ----------
        map_string : str
            Mapping formatted as ``"<bit>:<flag>, <bit>:<flag>, ..."``, e.g.
            ``"64:1, 2048:2"``

        Returns
        -------
        dict
            Mapping from healsparse bit value (int) to output flag value (int)

        Raises
        ------
        ValueError
            If an entry is not of the form ``<int>:<int>``

        """
        bit_flag_map = {}
        for entry in map_string.split(","):
            entry = entry.strip()
            if not entry:
                continue
            if not re.fullmatch(r"\d+\s*:\s*\d+", entry):
                raise ValueError(
                    f"Invalid BIT_FLAG_MAP entry '{entry}'; expected "
                    + "'<bit>:<flag>'"
                )
            bit, flag = (int(part) for part in entry.split(":"))
            bit_flag_map[bit] = flag
        if not bit_flag_map:
            raise ValueError("BIT_FLAG_MAP is empty")
        return bit_flag_map

    def _set_image_coordinates(self):
        """Set Image Coordinates.

        Read the image header into a WCS and record the image shape, mirroring
        ``Mask._set_image_coordinates``.

        """
        img = file_io.FITSCatalogue(self._image_path, hdu_no=0)
        img.open()
        self._header = img.get_header()
        # get_data().shape is (n_y, n_x)
        self._img_shape = img.get_data().shape
        img.close()
        del img

        self._wcs = wcs.WCS(self._header)

    def _default_chunk_size(self):
        """Default Chunk Size.

        Rows per chunk such that one band holds at most
        :data:`DEFAULT_PIXEL_BUDGET` pixels.

        Returns
        -------
        int
            Number of rows per chunk (at least 1)

        """
        n_x = self._img_shape[1]
        return max(1, DEFAULT_PIXEL_BUDGET // n_x)

    def _map_bits_to_flags(self, bit_values, off_map):
        """Map Bits to Flags.

        Translate healsparse bit values to output flag values through the
        bit→flag mapping, OR-combining every matching bit, and assign the
        off-map flag to sentinel pixels.

        Parameters
        ----------
        bit_values : numpy.ndarray
            Healsparse values queried at the pixel centres
        off_map : numpy.ndarray
            Boolean mask, ``True`` where the pixel lies outside the footprint

        Returns
        -------
        numpy.ndarray
            Output flag values (``int16``), same shape as ``bit_values``

        """
        flags = np.zeros(bit_values.shape, dtype=np.int16)
        for bit, flag in self._bit_flag_map.items():
            flags[(bit_values & bit) != 0] |= np.int16(flag)
        if self._off_map_flag != 0:
            flags[off_map] = np.int16(self._off_map_flag)
        return flags

    def rasterize(self):
        """Rasterize.

        Evaluate the pixel grid to (RA, Dec) in row chunks, query the
        healsparse mask, and build the ``int16`` flag image.

        Returns
        -------
        numpy.ndarray
            The rasterized flag image, shape ``(n_y, n_x)``, dtype ``int16``

        """
        # Lazy import: healsparse is an optional heavy dependency, imported at
        # use rather than module load.
        import healsparse

        hmap = healsparse.HealSparseMap.read(self._mask_path)
        sentinel = hmap.sentinel

        # Mask products come in two flavours: integer bit-flag maps (per-band
        # bits, e.g. 64 = r) and boolean maps (True = masked, e.g. the
        # candide copy of the 2025 r-band mask). For a boolean map the only
        # meaningful bit is 1 (True); any other bit silently selects nothing
        # (``True & 64 == 0`` -> an all-clean flag image), so fail loudly.
        if hmap.dtype == np.bool_:
            bad_bits = [bit for bit in self._bit_flag_map if bit != 1]
            if bad_bits:
                raise ValueError(
                    f"Mask {self._mask_path} is a boolean healsparse map; "
                    + f"BIT_FLAG_MAP bits {bad_bits} would never match "
                    + "(use '1:<flag>' for boolean masks)"
                )

        n_y, n_x = self._img_shape
        chunk_size = self._chunk_size or self._default_chunk_size()

        flag_image = np.zeros((n_y, n_x), dtype=np.int16)
        # Column indices are shared across every row band.
        x = np.arange(n_x)

        for y0 in range(0, n_y, chunk_size):
            y1 = min(y0 + chunk_size, n_y)
            yy, xx = np.meshgrid(np.arange(y0, y1), x, indexing="ij")
            # WCS uses 0-based pixel coordinates here (origin=0).
            ra, dec = self._wcs.all_pix2world(xx.ravel(), yy.ravel(), 0)
            # Normalize RA into [0, 360) so the wrap at 0/360 is handled;
            # healsparse expects lon in that range with lonlat=True.
            ra = np.mod(ra, 360.0)

            bit_values = hmap.get_values_pos(ra, dec, lonlat=True)
            off_map = bit_values == sentinel

            band = self._map_bits_to_flags(bit_values, off_map)
            flag_image[y0:y1, :] = band.reshape(y1 - y0, n_x)

        return flag_image

    def _combine_external_flag(self, flag_image):
        """Combine External Flag.

        Sum an external instrument flag image into the rasterized mask,
        matching ``Mask._build_final_mask`` semantics (element-wise sum of the
        two integer flag images).

        Parameters
        ----------
        flag_image : numpy.ndarray
            The rasterized flag image

        Returns
        -------
        numpy.ndarray
            The combined flag image (``int16``)

        """
        external_flag = file_io.FITSCatalogue(
            self._path_external_flag,
            hdu_no=self._hdu,
        )
        external_flag.open()
        ext_flag = external_flag.get_data()[:, :]
        external_flag.close()

        return (flag_image + ext_flag).astype(np.int16, copy=False)

    def _output_path(self):
        """Output Path.

        Full path of the output flag file, matching the ``mask`` module naming
        (``<PREFIX>_flag<num>.fits``).

        Returns
        -------
        str
            Output file path

        """
        name = (
            f"{self._img_prefix}{self._outname_base}"
            + f"{self._img_number}.fits"
        )
        return f"{self._output_dir}/{name}"

    def make_mask(self):
        """Make Mask.

        Rasterize the healsparse mask, optionally combine the external
        instrument flag, and write the flag file carrying the image WCS.

        Returns
        -------
        str
            Path to the written flag file

        """
        flag_image = self.rasterize()

        if self._path_external_flag is not None:
            flag_image = self._combine_external_flag(flag_image)

        output_path = self._output_path()
        out = file_io.FITSCatalogue(
            output_path,
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            hdu_no=0,
        )
        out.save_as_fits(
            data=flag_image,
            image=True,
            image_header=self._wcs.to_header(),
        )

        self._w_log.info(
            f"Wrote healsparse-derived flag file {output_path}"
        )

        return output_path
