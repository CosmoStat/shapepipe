"""MASK QUERY.

Class to flag SExtractor detections against external healsparse masks.

:Author: Claude Fable 5, for PR #847

"""

import numpy as np

from shapepipe.pipeline import file_io
from shapepipe.utilities import mask_query as mask_query_util


class MaskQuery(object):
    """Mask Query.

    Query external healsparse masks at every detection of a SExtractor
    catalogue and write the result as a single ``FLAG_EXT`` column into a copy
    of that catalogue.

    Parameters
    ----------
    sexcat_path : str
        Path to the input SExtractor catalogue
    output_path : str
        Path to the output catalogue
    mask_paths : list
        Paths to the healsparse maps to query
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
        ori_cat = file_io.FITSCatalogue(
            self._sexcat_path,
            SEx_catalogue=True,
        )
        ori_cat.open()
        data = ori_cat.get_data()
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
            "FLAG_EXT", flag, new_cat=True, new_cat_inst=new_cat
        )
        ori_cat.close()

        return int(np.count_nonzero(flag))
