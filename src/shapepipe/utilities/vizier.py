"""VIZIER QUERY UTILITY.

Consolidated Vizier query helper used by the mask module (to fetch the
reference star catalogue) and by ``scripts/python/create_star_cat.py``.
Retries over a list of mirror servers at progressively longer timeouts.

:Author: Martin Kilbinger

"""

import random
import time

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier


VIZIER_SERVERS = [
    "vizier.cds.unistra.fr",
    "vizier.cfa.harvard.edu",
    "vizier.iucaa.in",
]

VIZIER_TIMEOUTS = [10, 20, 40]


def query_vizier(ra, dec, radius_arcmin, cat_id):
    """Query a Vizier catalogue with retries over timeouts and mirror servers.

    Parameters
    ----------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    radius_arcmin : float
        Cone-search radius in arcminutes.
    cat_id : str
        Vizier catalogue identifier (e.g. ``"I/305/out"`` for GSC 2.3).

    Returns
    -------
    astropy.table.Table
        First result table returned by Vizier.

    Raises
    ------
    IndexError
        If all server/timeout combinations return an empty result.

    """
    # Empirically, single-precision input positions can cause Vizier to return
    # empty lists for some exposures; force double precision.
    p = np.array([ra, dec], dtype="double")
    coord = SkyCoord(ra=p[0] * u.deg, dec=p[1] * u.deg, frame="icrs")

    # Stagger concurrent queries to avoid hammering a single mirror.
    time.sleep(random.uniform(0, 5))

    for attempt, timeout in enumerate(VIZIER_TIMEOUTS):
        for server in VIZIER_SERVERS:
            v = Vizier(
                row_limit=-1, timeout=timeout, vizier_server=server
            )
            # cache=False: astroquery otherwise pickles every HTTP response into
            # $HOME/.astropy/cache/astroquery/Vizier, ~2 MB per query. The
            # workflow already caches the RESULT as a FITS catalogue on scratch
            # and skips the query when it hits, so the pickle is pure duplicate —
            # and at campaign scale (~25k exposures) it is ~50 GB against a
            # 50 GB home quota. Home is for source and config, not for a second
            # copy of the survey.
            result = v.query_region(
                coord, radius=radius_arcmin * u.arcmin, catalog=cat_id,
                cache=False,
            )
            if len(result) > 0:
                print(
                    f"Vizier query successful "
                    f"(server={server}, timeout={timeout}s)"
                )
                return result[0]
            print(
                f"Vizier returned empty list at {coord}, "
                f"{radius_arcmin:.2f} arcmin, "
                f"server={server}, timeout={timeout}s"
            )
        if attempt < len(VIZIER_TIMEOUTS) - 1:
            wait = 10 * 2**attempt
            print(
                f"All servers failed, retrying in {wait}s "
                f"(attempt {attempt + 1}/{len(VIZIER_TIMEOUTS)}, "
                f"next timeout={VIZIER_TIMEOUTS[attempt + 1]}s)"
            )
            time.sleep(wait)

    raise IndexError(
        f"Vizier astroquery returned empty list at {coord}, "
        f"radius={radius_arcmin} arcmin, catalog={cat_id}"
    )
