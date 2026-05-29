---
name: sqlitedict pickle-by-default is a known smell
status: open
tags:
    - shapepipe
    - security
    - future
    - refactor
created-at: 2026-05-28T09:37:09.591626506+02:00
outcome: Pickle-everywhere in SqliteDict storage (~10 modules) is dismissable today (trusted compute, no external write path) but a real smell. Revisit if shapepipe ever exposes a write path to untrusted .sqlite files.
---

## The vulnerability

`GHSA-g4r7-86gm-pgqc` (high, CVE-2024-35515) — sqlitedict uses pickle as the default value encoding. An attacker who can write `.sqlite` files that shapepipe later reads gets arbitrary code execution on `dict[key]` access, via crafted pickle payloads.

`first_patched_version: null`. The maintainer considers pickle-by-default intended behavior; the [README's serialization section](https://github.com/piskvorky/sqlitedict?tab=readme-ov-file#serialization) tells users to pass `encode=`/`decode=` if they care. There will be no upstream fix.

## How shapepipe uses it

`SqliteDict(path)` with default (pickle) encoding, across ~10 modules:

- `make_cat`, `psfex_interp`, `mccd_interp`, `ngmix`, `vignetmaker`, `merge_headers`, `sextractor`, `python_example_runner`
- Stores numpy arrays, astropy WCS objects, astropy Header objects, and dicts of those — intermediate pipeline state passed between runner stages

None of the call sites override `encode=` / `decode=`. All use the vulnerable default.

## Why dismissable today

Threat model requires attacker-controlled `.sqlite` files reaching a pipeline run. In normal shapepipe operation those files are produced and consumed by the same pipeline run, on the same trusted compute node (candide / canfar / cineca / feynman). No external write path. An attacker who can write to a shapepipe scratch dir already has filesystem access to a trusted compute node — at which point poisoning a vignet cache is the least of our worries.

Dismissed in Dependabot UI on 2026-05-28 with `tolerable_risk` reason, pointing at this fiber.

## Why the encode/decode workaround isn't quite worth it

Considered three paths to actually close the vector:

| Option | Effort | Catches | Catch |
|---|---|---|---|
| HMAC-pickle (`encode=hmac_pickle`, `decode=hmac_verify_then_unpickle`) | ~4h | File-tampering by attackers without the HMAC key | Only useful if the attacker can't read the key — but they're co-located on the compute node, so they can. Closes nothing actually relevant. |
| msgpack + msgpack-numpy + custom astropy hooks | 1–2 days | All RCE via deserialization | WCS/Header round-trip is fragile; hooks often end up calling pickle internally for unknown astropy types. Brittle and hard to verify. |
| Replace with HDF5/parquet across the pipeline | Days–week | The whole class of issues | Real refactor of ~10 modules; on-disk format change breaks existing scratch files mid-run. The right architectural answer but not a CVE-triage move. |

The first two don't actually move the security frontier given the threat model. The third is real work that should land for architectural reasons (typed, schema-validated intermediates), not because of this CVE.

## When to revisit

- shapepipe gains a web-service / API surface where users upload data
- Shared workspaces where multiple non-coordinated users write to the same scratch dirs
- Any deployment where the trust boundary between pipeline-input and pipeline-output narrows
- The HDF5/parquet intermediate-state refactor lands for other reasons (consolidating typing / schemas) — then close this fiber as superseded

## Adjacent: dismissed-alert hygiene

This is the first `tolerable_risk` dismissal on shapepipe. If more accumulate, consider grouping the analysis into a single `shapepipe/security-posture` reference fiber that links out to each dismissed alert's reasoning. One fiber per dismissed CVE is fine for now (low volume), pyramidalizes to a posture doc later.
