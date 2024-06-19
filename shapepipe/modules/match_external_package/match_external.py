"""MATCH EXTERNAL.

This module matches an external catalogue to a ShapePipe (SExtractor)
catalogue.

:Authors: Martin Kilbinger, Xavier Jimenez

"""

import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord, match_coordinates_sky

from shapepipe.pipeline import file_io


def get_cat(path):
    """Get Catalogue.

    Open a FITS catalogue.

    Parameters
    ----------
    path : str
        Path to catalogue

    Returns
    -------
    file_io.FITSCatalogue
        Open FITS catalogue object

    """
    cat = file_io.FITSCatalogue(path)
    cat.open()

    return cat


def get_data(path, hdu_no):
    """Get Data.

    Extract data from a given catalogue HDU.

    Parameters
    ----------
    path : str
        Path to catalogue
    hdu_no : int
        HDU number

    Returns
    -------
    tuple
        Data, column names and extension names

    """
    cat = get_cat(path)
    data = cat.get_data(hdu_no)
    col_names = cat.get_col_names(hdu_no=hdu_no)
    ext_names = cat.get_ext_name()
    cat.close()

    return data, col_names, ext_names


def get_ra_dec(data, col_ra, col_dec):
    """Get RA and Dec.

    Get RA and Dec from input data array.

    Parameters
    ----------
    data : numpy.ndarray
        Input data array
    col_ra : int
        Column number for RA
    col_dec : int
        Column number for Dec

    Returns
    -------
    tuple
        RA and Dec values

    """
    ra = data[col_ra]
    dec = data[col_dec]

    return ra, dec


class MatchCats(object):
    """Match Catalogues.

    Parameters
    ----------
    input_file_list : list
        List of input catalogue paths to be pasted
    output_path : str
        Output file path of pasted catalogue
    w_log : logging.Logger
        Logging instance
    tolerance : astropy.units.quantity.Quantity
        Tolerance in arcsec
    col_match : list
        (Internal data) column name(s) to copy into matched output catalogue
    hdu_no : int
        (Internal) catalogue HDU number
    mode : str
        Run mode, ``CLASSIC`` or ``MULTI-EPOCH``
    external_cat_path : str
        External catalogue path
    external_col_match : list
        External data column name(s) for matching
    external_col_copy : list
        Column name(s) to copy into matched output catalogue
    external_hdu_no : int, optional
        External catalogue hdu number, default is ``1``
    mark_non_matched : float, optional
        If not ``None``, output not only matched but all objects, and mark
        non-matched objects with this value
    output_distance : bool, optional
        Output distance between matches if ``True``, default is ``False``

    """

    def __init__(
        self,
        input_file_list,
        output_path,
        w_log,
        tolerance,
        col_match,
        hdu_no,
        mode,
        external_cat_path,
        external_col_match,
        external_col_copy,
        external_hdu_no=1,
        mark_non_matched=None,
        output_distance=False,
    ):

        self._input_file_list = input_file_list
        self._output_path = output_path
        self._w_log = w_log

        self._tolerance = tolerance * units.arcsec

        self._col_match = col_match
        self._hdu_no = hdu_no
        self._mode = mode

        self._external_cat_path = external_cat_path
        self._external_col_match = external_col_match
        self._external_col_copy = external_col_copy
        self._external_hdu_no = external_hdu_no

        self._mark_non_matched = mark_non_matched
        self._output_distance = output_distance

    def process(self):
        """Process.

        Process catalogues.

        """
        # Load external and internal data
        external_data, dummy1, dummy2 = get_data(
            self._external_cat_path,
            self._external_hdu_no,
        )
        external_ra, external_dec = get_ra_dec(
            external_data,
            self._external_col_match[0],
            self._external_col_match[1],
        )
        external_coord = SkyCoord(ra=external_ra, dec=external_dec, unit="deg")

        data, col_names, ext_names = get_data(
            self._input_file_list[0],
            self._hdu_no,
        )
        ra, dec = get_ra_dec(data, self._col_match[0], self._col_match[1])
        coord = SkyCoord(ra=ra, dec=dec, unit="deg")

        # Match objects in external cat to internal cat. indices=indices to
        # external object for each object in internal cat e.g.
        # external_coord[indices[0]] is the match for coord[0].
        indices, d2d, d3d = match_coordinates_sky(
            coord,
            external_coord,
            nthneighbor=1,
        )

        # Find close neighbours, indices_close is True for all close matches
        indices_close = d2d < self._tolerance

        if not any(indices_close):
            self._w_log.info(
                f"No match for {self._input_file_list[0]} with distance < "
                + f"{self._tolerance} arcsec found, no output created."
            )

        else:
            # Get indices in internal and external catalogues of pair-wise
            # matches
            w = np.array(
                [
                    (idx, ide)
                    for (idx, ide) in enumerate(indices)
                    if indices_close[idx]
                ]
            )
            id_sub = w[:, 0]
            id_ext_sub = w[:, 1]
            id_all = np.arange(len(indices))

            if self._mark_non_matched:
                # Output all objects
                id_data = id_all
                id_ext = indices
            else:
                # Output only matched objects
                id_data = id_sub
                id_ext = id_ext_sub

            self._w_log.info(
                f"{len(id_sub)} objects matched out of {len(indices)}."
            )

            # Copy matched objects from internal catalogue to output data
            matched = {}
            for col in col_names:
                matched[col] = data[col][id_data]

            # Copy columns from external catalogue to output data
            for col in self._external_col_copy:
                matched[col] = external_data[col][id_ext]
                if self._mark_non_matched:
                    for idx, i_ext in enumerate(indices):
                        if not indices_close[idx]:
                            matched[col][idx] = self._mark_non_matched

            # Output distance if desired
            if self._output_distance:
                # Output distance in arcsec
                matched["distance"] = d2d[id_data].to("arcsec").value

            # Write FITS file
            out_cat = file_io.FITSCatalogue(
                self._output_path,
                SEx_catalogue=False,
                open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            )
            out_cat.save_as_fits(
                data=matched,
                ext_name="MATCHED",
            )

            # Write all extensions if in multi-epoch mode
            if self._mode == "MULTI-EPOCH":
                hdu_me_list = [
                    idx for idx, name in enumerate(ext_names) if "EPOCH" in name
                ]
                for hdu_me in hdu_me_list:
                    data_me, col_names_me, dummy = get_data(
                        self._input_file_list[0],
                        hdu_me,
                    )
                    matched_me = {}
                    for col_me in col_names_me:
                        matched_me[col_me] = data_me[col_me][id_data]
                    out_cat.save_as_fits(
                        data=matched_me,
                        ext_name=ext_names[hdu_me],
                    )
