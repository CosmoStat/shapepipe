"""STAR COLOUR.

Attach an external colour to each PSF star and drop stars an external purity
verdict rejects, between star selection and the PSF fit.

The external catalogue is a FITS table with one row per object: world
coordinates, a colour (g−i), and a boolean rejection flag. Each star in the
training and reserve splits is matched to its nearest catalogue object within
the match radius. Training stars whose match is flagged are dropped; the
reserve is kept whole so that it stays comparable with an uncoloured run. Every
star gets a ``GI_COLOR`` column for the PSF model's context: its matched colour,
or, for a star without one, the median over the CCD's kept training stars.

"""

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree


def _xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


class StarColour:
    """Match PSF stars to an external colour catalogue.

    @sc [decision:star_selection_psf.psf_colour_purification]

    Parameters
    ----------
    catalogue_path : str
        FITS table with ``RA``, ``DEC`` (deg), the colour column and the
        rejection column
    colour_col, reject_col : str
        Column names in the catalogue
    radius_arcsec : float
        Match radius
    """

    def __init__(self, catalogue_path, colour_col="GI", reject_col="REMOVE", radius_arcsec=0.5):
        cat = fits.getdata(catalogue_path, 1)
        self._colour = np.asarray(cat[colour_col], dtype=float)
        self._reject = np.asarray(cat[reject_col], dtype=bool)
        self._tree = cKDTree(_xyz(cat["RA"], cat["DEC"]))
        self._radius = np.radians(radius_arcsec / 3600)

    def match(self, ra, dec):
        """Return colour (NaN if unmatched or colourless) and rejection flag."""
        d, i = self._tree.query(_xyz(ra, dec), distance_upper_bound=self._radius)
        ok = np.isfinite(d)
        j = np.where(ok, i, 0)
        return np.where(ok, self._colour[j], np.nan), ok & self._reject[j]

    @staticmethod
    def _write(src, dst, keep, colour):
        h = fits.open(src)
        data = h[2].data[keep]
        cols = [fits.Column(name=c.name, format=c.format, unit=c.unit, dim=c.dim, array=np.asarray(data[c.name]))
                for c in h[2].columns]
        cols.append(fits.Column(name="GI_COLOR", format="1E", array=colour[keep].astype(np.float32)))
        fits.HDUList([h[0], h[1], fits.BinTableHDU.from_columns(cols, name="LDAC_OBJECTS")]).writeto(dst, overwrite=True)

    def process(self, train_path, reserve_path, train_out, reserve_out, w_log):
        """Write the purified, coloured training split and the coloured reserve."""
        tr = fits.getdata(train_path, 2)
        c_tr, rej = self.match(tr["XWIN_WORLD"], tr["YWIN_WORLD"])
        have = np.isfinite(c_tr) & ~rej
        fill = float(np.median(c_tr[have])) if have.any() else 1.3
        self._write(train_path, train_out, ~rej, np.where(np.isfinite(c_tr), c_tr, fill))

        rs = fits.getdata(reserve_path, 2)
        c_rs, _ = self.match(rs["XWIN_WORLD"], rs["YWIN_WORLD"])
        self._write(reserve_path, reserve_out, np.ones(len(rs), bool), np.where(np.isfinite(c_rs), c_rs, fill))
        w_log.info(f"star_colour: {rej.sum()} of {len(tr)} training stars rejected; "
                   f"{(~np.isfinite(c_tr)).sum()} without colour given {fill:.3f}")
