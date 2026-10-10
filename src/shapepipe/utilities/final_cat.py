"""FINAL CATALOGUE.

Read a per-tile final catalogue, the output of ``make_cat_runner``, in either
of the formats it exists in:

* HDF5, as ``make_cat.write_final_cat`` writes it: one dataset per column, in
  column order, a vector column as a 2-D dataset with one row per object;
* FITS, a binary table (``RESULTS``, HDU 1).

The format is decided per file from its content (the HDF5 signature), not its
name, so callers can read a directory containing both formats.

"""

import h5py
import numpy as np


def is_hdf5(path):
    """Is HDF5.

    Parameters
    ----------
    path : str
        Catalogue path

    Returns
    -------
    bool
        ``True`` if the file carries the HDF5 signature, ``False`` otherwise
        (taken to be FITS)

    """
    return h5py.is_hdf5(path)


def read_final_cat(path, columns=None, hdu=1):
    """Read Final Catalogue.

    Read a per-tile final catalogue, HDF5 or FITS, into a structured array.
    Every field is in native byte order and a vector column is a sub-array
    field, so the same catalogue reads to the same array from either format.

    Parameters
    ----------
    path : str
        Catalogue path
    columns : list of str, optional
        Columns to read, in the order wanted; ``None`` or empty reads every
        column, in file order
    hdu : int, optional
        HDU of the table in a FITS catalogue; default is ``1``

    Returns
    -------
    numpy.ndarray
        Structured array, one field per column

    Raises
    ------
    KeyError
        If the catalogue lacks a requested column; the message names the file
        and every missing column

    """
    if is_hdf5(path):
        with h5py.File(path, "r") as cat:
            present = list(cat.keys())
            wanted = list(columns) if columns else present
            _check_columns(path, wanted, present)
            arrays = {col: cat[col][()] for col in wanted}
    else:
        from astropy.io import fits

        with fits.open(path, memmap=False) as hdu_list:
            data = hdu_list[hdu].data
            present = list(data.columns.names)
            wanted = list(columns) if columns else present
            _check_columns(path, wanted, present)
            arrays = {col: np.asarray(data[col]) for col in wanted}

    return _to_structured(arrays)


def _check_columns(path, wanted, present):
    """Raise ``KeyError`` naming every column of ``wanted`` not in ``present``."""
    present = set(present)
    missing = [col for col in wanted if col not in present]
    if missing:
        raise KeyError(
            f"{path}: missing {len(missing)} of the {len(wanted)} requested "
            f"column(s): {' '.join(missing)}"
        )


def _to_structured(arrays):
    """Assemble column arrays, one row per object, into a structured array."""
    dtype = np.dtype(
        [
            (col, a.dtype.newbyteorder("="), a.shape[1:])
            for col, a in arrays.items()
        ]
    )
    lengths = {len(a) for a in arrays.values()}
    if len(lengths) > 1:
        raise ValueError(f"columns differ in length: {sorted(lengths)}")
    out = np.empty(lengths.pop() if lengths else 0, dtype=dtype)
    for col, a in arrays.items():
        out[col] = a
    return out
