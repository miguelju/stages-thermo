# Notebooks — the learning path

One executable notebook per milestone, following the vle notebook conventions
(setup cell, context, worked example, ≥2 exercises with hidden solutions, pinned
assertion cells). The rigorous-solver notebooks show the machinery — residuals,
Jacobian sparsity, convergence history — not just the answer ("open the black
box," PLAN §1).

Planned (see `ROADMAP.md` for status):

| Notebook | Milestone |
|---|---|
| [`01-mccabe-thiele.ipynb`](01-mccabe-thiele.ipynb) | **M1 — live** — benzene–toluene (PR), methanol–water (van Laar, ties to vle Ch. IV); 3 exercises |
| [`02-ponchon-savarit.ipynb`](02-ponchon-savarit.ipynb) | **M2 — live** — benzene–toluene consistency check, methanol–water CMO drift, acetone–water CMO design failure (curved operating lines, R_min underestimated), ammonia–water two-route showcase on a Pátek–Klomfar-digitized reference chart; 2 exercises |
| `03-shortcut-design.ipynb` | M3 — FUG depropanizer design |
| `04-mesh-and-bubble-point.ipynb` | M6 — "how a column is actually solved," part 1 |
| `05-absorbers-sum-rates.ipynb` | M7 |
| `06-naphtali-sandholm.ipynb` | M8 — Jacobian sparsity, quadratic convergence, part 2 |
| `07-when-newton-fails.ipynb` | M9 — homotopy + pseudo-transient |
| `08-numpy-for-columns.ipynb` | M10 — batch API, N-vs-R maps |
| `09-inside-out.ipynb` | M11 (stretch) |

Notebooks drive **both** packages: `vle.System` for thermodynamics exploration
and `stages` for columns.

Execute before committing (per CLAUDE.md notebook conventions):

```sh
~/miniconda3/envs/stages/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/01-mccabe-thiele.ipynb
```
