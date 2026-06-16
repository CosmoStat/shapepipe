---
id: 01KTT34NW92KPHH9GDDXAGZZ89
name: 'Fresh-pass review: PR #741 (post-integration)'
tags:
    - shapepipe
    - ngmix
    - review
created-at: 2026-06-11T03:00:58.633703716+02:00
outcome: 'Fresh full-diff review of #741 at b2dcd793 (empirical, against ngmix 2.4.0): scientific core verified sound (seed-reproducibility holds; injected shear recovered to 0.2% with response correction; r50/ivar-weights math checks out) but two empirically demonstrated blockers found in the v2.0 rewrite''s robustness shell — TILE_ID metadata key truncates the post-process CCD scan (~80% of epochs silently lost to multi-epoch bookkeeping) and failed galaxy fits crash the whole tile at save (KeyError, v1''s per-object BootGalFailure contract lost). Plus: one failed PSF epoch kills the object despite ignore_failed_psf; mcal_flags hardcoded 0 (dead quality column); cfis_simu configs broken against the positional-WCS contract. Science note: with psf=''fitgauss'' the reconvolved metacal PSF is round by construction, so g1/g2_psfo are identically zero — PSF-leakage regressions on them are degenerate. Fixes + review part 3 posted to the PR same night.'
---

Second-look review requested by Cail after the constitution-sweep integration round, run by a fresh agent with no investment in the branch. Full findings live in PR #741's "Review — part 3" comment; this fiber is the pointer plus what matters beyond the PR.

## Why the blockers hid

Both blockers survive a green unit-test suite and a one-tile smoke run with easy objects: the TILE_ID truncation only bites multi-CCD bookkeeping breadth (the first 8 CCDs still work), and the compile_results KeyError needs a hard-to-fit object (easy smoke objects all converge). Lesson for the test suite: structural tests catch import/config rot; these needed *adversarial empirical* tests — forced-failure objects, metadata-polluted sqlite fixtures. The regression tests added with the fixes are exactly that shape.

## Science threads left deliberately open

- **fixnoise homoscedasticity** (#604 thread): metacal noise-symmetrization currently draws at median(bkg_rms) while weights are per-pixel — second-order within a 51px stamp, one-line to change, Martin's call.
- **g1/g2_psfo ≡ 0 under psf='fitgauss'** (*_psfo thread): keep/rename/drop decision affects any downstream PSF-leakage regression; relevant to the tutorial.
