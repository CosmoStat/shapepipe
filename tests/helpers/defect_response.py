"""Full-matrix metacal recovery for fixed detector defects and neighbours."""

import galsim
import numpy as np

from shapepipe.modules.ngmix_package import ngmix as ngm
from shapepipe.testing.simulate import make_data
from tests.helpers.metacal_sim import build_stamp


N = 51
PIXEL_SCALE = 0.1857
FLUX = 1000.0


def neighbour_scene(offset, shear, psf, psf_shear, hlr=0.5, flux_ratio=1.0,
                    seg_radius=3):
    """A neighbour's image and the stamp's segmentation map.

    The neighbour is an exponential galaxy whose centre lies ``offset``
    (row, col) pixels from the stamp centre, sheared with the scene (its
    profile and its position) and convolved with the detector-frame PSF.
    The segmentation map labels a disk of ``seg_radius`` px around the
    stamp centre 1 (the target) and one around the unsheared neighbour
    position 2.
    """
    row, col = offset
    matrix = galsim.Shear(g1=shear[0], g2=shear[1]).getMatrix()
    dx, dy = matrix @ np.array([col, row]) * PIXEL_SCALE
    image = galsim.Convolve(
        galsim.Moffat(beta=2.5, fwhm=psf).shear(
            g1=psf_shear[0], g2=psf_shear[1]),
        galsim.Exponential(half_light_radius=hlr, flux=FLUX * flux_ratio)
        .shear(g1=shear[0], g2=shear[1]),
    ).shift(dx, dy).drawImage(nx=N, ny=N, scale=PIXEL_SCALE).array
    rows, cols = np.indices((N, N)) - N // 2
    seg = np.zeros((N, N), dtype=np.int32)
    seg[rows ** 2 + cols ** 2 <= seg_radius ** 2] = 1
    seg[(rows - row) ** 2 + (cols - col) ** 2 <= seg_radius ** 2] = 2
    return image, seg


def defect_response(bad, hlr=0.5, psf=0.7, seeds=range(4),
                    psf_shear=(0.0, 0.0), options=None, known_rms=True,
                    defect_value=1e3, neighbour=None):
    """Recover c and M = inverse(mean R) A - I with paired seed errors.

    Null pairs rotate the pixels by 90 degrees while pre-rotating the PSF
    ellipticity oppositely, so the final PSF stays fixed in detector
    coordinates.
    This does NOT average away elliptical-PSF leakage.

    The flagged pixels ``bad`` hold ``defect_value`` (far above the galaxy's
    peak) rather than sky, as a hot column or bleed would, so any defect value
    that reaches metacal shows up as a bias.

    ``neighbour``, if given, holds the keyword arguments of
    :func:`neighbour_scene` but ``shear``, ``psf`` and ``psf_shear``: a
    neighbour fixed in detector coordinates like the defect. Its
    segmentation footprint is the stamp's neighbour mask under
    ``blend_handling="noisefill"``, and the segmentation map rides the stamp
    for ``"uberseg"`` (pass ``object_number=1`` in ``options``).

    The result also holds ``s2n``, the mean over arms of the no-shear fit's
    S/N, the price of the pixels a treatment drops.
    """
    gamma = 0.02
    options = {} if options is None else options
    samples, s2n = [], []
    for seed in seeds:
        def arm(axis, sign, rotation=0):
            shear = [0.0, 0.0]
            if axis >= 0:
                shear[axis] = sign * gamma
            ps = tuple(v * (-1 if rotation else 1) for v in psf_shear)
            data = list(make_data(
                rng=np.random.RandomState(seed + 100), shear=shear,
                psf_shear=ps, noise=1e-4, n_epochs=1, img_size=51,
                gal_hlr=hlr, psf_fwhm=psf, return_centers=True,
            ))
            centre = data.pop()[0]
            offset = np.array([centre.y - 26, centre.x - 26])
            if rotation:
                data[0] = [np.rot90(a).copy() for a in data[0]]
                data[1] = [np.rot90(a).copy() for a in data[1]]
                offset = np.array([-offset[1], offset[0]])
            seg = None
            if neighbour is not None:
                light, seg = neighbour_scene(
                    shear=shear, psf=psf, psf_shear=psf_shear, **neighbour)
                data[0] = [a + light for a in data[0]]
            data[0] = [np.where(bad, defect_value, a) for a in data[0]]
            data[4] = [bad.astype(np.int32)]
            stamp = build_stamp(data)
            stamp.offsets = [offset]
            if seg is not None:
                stamp.segs = [seg]
                stamp.neighbours = [seg == 2]
            if known_rms:
                stamp.bkg_rms = [np.full((51, 51), 1e-4)]
            rng = np.random.RandomState(seed)
            result, _, _ = ngm.do_ngmix_metacal(
                stamp, ngm.get_prior(0.1857, rng), 1.0, rng,
                centroid_source="wcs", **options,
            )
            assert all(result[t]["flags"] == 0 for t in ngm.METACAL_TYPES)
            s2n.append(result["noshear"]["s2n"])
            e = np.asarray(result["noshear"]["g"])
            response = np.column_stack([
                (np.asarray(result[p]["g"]) - result[m]["g"]) / 0.02
                for p, m in (("1p", "1m"), ("2p", "2m"))
            ])
            assert np.all(np.isfinite(e)) and np.all(np.isfinite(response))
            return e, response
        null = [arm(-1, 0, k) for k in (0, 1)]
        e0 = np.mean([a[0] for a in null], axis=0)
        r0 = np.mean([a[1] for a in null], axis=0)
        derivatives, responses = [], []
        for axis in (0, 1):
            plus, minus = arm(axis, 1), arm(axis, -1)
            derivatives.append((plus[0] - minus[0]) / (2 * gamma))
            responses.extend([plus[1], minus[1]])
        samples.append((e0, r0, np.column_stack(derivatives),
                        np.mean(responses, axis=0)))
    e, r, a, rm = [np.array([s[i] for s in samples]) for i in range(4)]
    assert np.linalg.svd(rm.mean(axis=0), compute_uv=False).min() > 0.1
    c = np.linalg.solve(r.mean(axis=0), e.mean(axis=0))
    matrix = np.linalg.solve(rm.mean(axis=0), a.mean(axis=0)) - np.eye(2)
    draw = np.random.RandomState(91).randint(len(e), size=(1000, len(e)))
    cb = np.linalg.solve(r[draw].mean(axis=1),
                         e[draw].mean(axis=1)[..., None])[..., 0]
    mb = np.linalg.solve(rm[draw].mean(axis=1), a[draw].mean(axis=1))
    mb -= np.eye(2)
    return dict(c=c.tolist(), m=np.diag(matrix).tolist(),
                matrix=matrix.tolist(),
                c_err=cb.std(axis=0).tolist(),
                m_err=np.diag(mb.std(axis=0)).tolist(),
                response=r.mean(axis=0).tolist(),
                s2n=float(np.mean(s2n)),
                seed_groups=[dict(e=s[0].tolist(), R=s[1].tolist(),
                                  A=s[2].tolist(), Rm=s[3].tolist())
                             for s in samples])
