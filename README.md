# RevascLogic — The Completeness Paradox

A single-file, **fully offline** dashboard that pools odds ratios for
**culprit-only vs complete revascularization** in acute coronary syndromes,
using a Hartung–Knapp–Sidik–Jonkman (HKSJ) random-effects model. It contrasts an
efficacy signal (recurrent MI / MACE) against a hard-safety signal
(death / renal), shows a forest plot, a 100-patient competing-risks waffle, and a
Monte-Carlo net-clinical-benefit distribution, with preset scenarios
(Stable-only, JACC-mixed, Shock-only, All-comers).

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN.

## Layout

```
index.html   single-file UI (loads engine.js)
engine.js    pure statistical core — runs in Node and the browser
tests.js     Node test harness, 28 assertions
LICENSE      Apache-2.0
```

## Statistical core (`engine.js`)

| Function | What it does |
|---|---|
| `runHKSJ(data, kEff, kSe)` | HKSJ random-effects pooling of log odds ratios on the supplied effect/SE keys: DerSimonian–Laird τ² and I² via `(Q−(k−1))/C`, HKSJ variance floored at the classical RE variance, and a `t_{k−1}` critical value. Returns `{ mu, se_pool, tVal, tau2, I2, trials }` or `null` on empty input |
| `tInv975(df)` | exact two-sided 0.975 Student-t critical value (table for df 1..30, Cornish–Fisher for larger df) |

Pooling is done on the **log** scale and exponentiated for display. The
stochastic competing-risk / net-clinical-benefit Monte-Carlo projection stays in
the page (`computeResults` in `index.html`) because it depends on `Math.random`;
all deterministic pooling is in `engine.js` as the single source of truth.

## Fixes applied during revival (2026-06-05)

1. **Offline:** removed the Google Fonts `<link>`; the app now loads no external
   resource (system fonts fall back).
2. **Single source of truth:** extracted the HKSJ core into a pure `engine.js`.
   The former web-worker (a Blob built from inline text) was removed and the
   pooling now runs synchronously; the inline duplicate of `runHKSJ` was deleted
   and the page loads `engine.js`.
3. **HKSJ variance floor (correctness bug).** The original computed the HK
   variance as `q_term / ((k−1)·ΣwRE)` with **no floor**, so when `Q < k−1`
   (HK multiplier < 1) the confidence interval narrowed *below* the classical
   random-effects interval — anti-conservative, and for homogeneous studies it
   collapsed the CI toward a point. It is now floored at the RE variance:
   `var_HKSJ = max(1, qStat)·(1/ΣwRE)`, the standard Hartung–Knapp form. This is
   the well-documented HKSJ floor and is verified by the two-identical-study test.
4. **Proper t critical value (correctness bug).** The original hardcoded crude
   approximations (`k<10 → 2.3`, etc.). It now uses an exact `qt(0.975, k−1)`
   table (e.g. k=2 → 12.706, k=3 → 4.303, k=5 → 2.776).

The DerSimonian–Laird τ²/I²/Q math was verified correct and left unchanged.

## Tests

```
node tests.js
# 28 passed, 0 failed
```

Every expected value is hand-derived independently in `tests.js` (not produced by
the engine). Checks include: the `tInv975` table, empty and `null` guards, a
single-study passthrough, a two-identical-study case exercising the HKSJ floor
(τ²=0, I²=0, se_pool=√0.02), a hand-worked heterogeneous two-study pooling
(τ²≈0.0855, pooled OR≈0.653, I²≈77.4%, se_pool≈0.2328), and per-trial CI bounds.

## Caveats

DerSimonian–Laird is known to under-estimate τ² for small *k* (REML/Paule–Mandel
are preferred for k<10); the dashboard preserves the original method for
continuity and reports τ² and I² alongside every estimate. The trial database is
an illustrative reconstruction and the Monte-Carlo projection is a simplified
competing-risks abstraction, so results are hypothesis-generating, not a clinical
decision rule. Apache-2.0 licensed.
