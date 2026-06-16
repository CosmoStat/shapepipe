---
id: 01KTK3WGKFCGD2VQAQSABC3TBJ
name: np.str0 fix (NumPy 2 aftershock)
status: closed
tags:
    - shapepipe
    - numpy
    - bugfix
created-at: 2026-06-08T09:59:18.639370578+02:00
closed-at: 2026-06-10T17:13:38.507495758+02:00
outcome: 'np.str0 removed from file_io.py on develop (verified 2026-06-10: only np.str_ remains); NumPy 2 string-column saves work. Regression guard idea (NAME column in roundtrip test) noted in body if ever needed.'
---

Fabian Hervas-Peters hit `module 'numpy' has no attribute 'str0'` running the
new container against `src/shapepipe/pipeline/file_io.py`. Root cause: the
container now ships **NumPy 2.4.4** (`pyproject` pins `numpy>=2.0`), and
`np.str0` — a deprecated alias for `np.str_` — was *removed* in NumPy 2.0.

The offending line was `_get_fits_col_type`'s string branch:
`elif type(col_data[0]) in [str, np.str_, np.str0]:`. The fix is to drop the
dead alias — `np.str_` already covers it, so nothing is lost.

**Why the test suite never caught it.** The list literal `[str, np.str_, np.str0]`
is only *evaluated* when that `elif` is actually reached, i.e. only when a column
whose first element is a string flows through `save_as_fits`. The existing
`test_fits_catalogue_table_roundtrips` saved only a `float64` and an `int16`
column — both short-circuit on earlier branches and never touch the string path.
So the `AttributeError` stayed dormant until Fabian saved a catalogue with a
string column. The regression guard is a one-line addition: a `NAME` string
column in the roundtrip test, which now exercises the `"A"` (FITS char) dispatch.

Interestingly, `ngmix_v2.0` had *already* removed `np.str0` independently, so
merging develop into it was a clean no-op on `file_io.py` and only carried over
the improved test.

General lesson: removed-NumPy-2 aliases (`np.str0`, `np.bool0`, `np.int0`,
`np.object0`, `np.unicode_`, …) hide in rarely-reached branches and survive
import — they only blow up at the call site. Worth a grep sweep when migrating a
codebase to NumPy 2. A sweep of `src/shapepipe/` on 2026-06-08 found `str0` was
the only remaining offender.
