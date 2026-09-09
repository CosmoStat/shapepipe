"""MASK QUERY.

Class to flag SExtractor detections against external healsparse masks.

:Author: Claude Fable 5, for PR #847

"""

import shutil

import numpy as np

from shapepipe.pipeline import file_io
from shapepipe.utilities import mask_query as mask_query_util


class MaskQuery(object):
    """Mask Query.

    Query external healsparse masks at every detection of a SExtractor
    catalogue and write the result as a single ``MASK_EXT`` column into a copy
    of that catalogue.

    With no maps configured this is a strict no-op: the input catalogue is
    copied through unchanged, with no ``MASK_EXT`` column, matching
    ``make_cat``'s ``MASK_EXT_PATHS`` contract (absent key, nothing happens).
    The copy is what keeps the module in the chain — ``setools`` reads this
    module's output, so producing no file would break the chain rather than
    disable the query.

    Parameters
    ----------
    sexcat_path : str
        Path to the input SExtractor catalogue
    output_path : str
        Path to the output catalogue
    mask_paths : list
        Paths to the healsparse maps to query; empty means no-op
    bits : int, optional
        Bit mask applied to integer maps; default ``None`` flags any nonzero
        value
    w_log : logging.Logger, optional
        Logging instance

    """

    def __init__(
        self,
        sexcat_path,
        output_path,
        mask_paths,
        bits=None,
        w_log=None,
    ):

        self._sexcat_path = sexcat_path
        self._output_path = output_path
        self._mask_paths = mask_paths
        self._bits = bits
        self._w_log = w_log

    def process(self):
        """Process.

        Query the masks and write the flagged catalogue.

        Returns
        -------
        int
            Number of flagged objects

        """
        if not self._mask_paths:
            # No maps configured: copy the catalogue through byte-for-byte.
            # A rewrite through FITSCatalogue would round-trip the LDAC HDUs
            # for no reason; this way an unconfigured run is provably
            # identical to its input.
            shutil.copyfile(self._sexcat_path, self._output_path)
            if self._w_log is not None:
                self._w_log.info(
                    "No MASK_PATHS configured; passing "
                    + f"{self._sexcat_path} through unchanged, no MASK_EXT"
                    + " column written"
                )
            return 0

        ori_cat = file_io.FITSCatalogue(
            self._sexcat_path,
            SEx_catalogue=True,
        )
        ori_cat.open()
        data = ori_cat.get_data()

        # A CCD SExtractor found nothing on is tolerated all along this chain
        # (setools' ~0.2% attrition, psfex_interp's floor=0 warn), so it must
        # not be an error here either. An empty LDAC table has no columns to
        # index, so read the positions only when there are rows, and still
        # publish an output file — a missing sexcat_ext would look to the file
        # handler like a crash rather than like an empty CCD.
        if len(data) == 0:
            ra = np.zeros(0)
            dec = np.zeros(0)
            if self._w_log is not None:
                self._w_log.info(
                    "No detections in "
                    + f"{self._sexcat_path}; writing an empty MASK_EXT column"
                )
        else:
            ra = np.copy(data["XWIN_WORLD"])
            dec = np.copy(data["YWIN_WORLD"])

        flag = mask_query_util.flag_positions(
            self._mask_paths,
            ra,
            dec,
            bits=self._bits,
            w_log=self._w_log,
        )
        # int32 is what the catalogue carries; the UNIONS bit table needs 12
        # bits, and no OR of it can overflow.
        flag = flag.astype(np.int32)

        new_cat = file_io.FITSCatalogue(
            self._output_path,
            SEx_catalogue=True,
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
        )
        ori_cat.add_col(
            "MASK_EXT", flag, new_cat=True, new_cat_inst=new_cat
        )
        ori_cat.close()

        return int(np.count_nonzero(flag))
