# E156-PROTOCOL — RevascLogic (The Completeness Paradox)

- **Project:** culpritcontroversy (GitHub repo `culpritcontroversy`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `culpritcotro.html` dump)
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: removed the Google Fonts CDN `<link>`; the app now
  loads no external resource (system fonts fall back).
- Extracted the statistical core into a pure `engine.js` (single source of
  truth); the former Blob web-worker was removed and pooling now runs
  synchronously, with the inline duplicate of `runHKSJ` deleted.
- Added `tests.js` (28 assertions, all passing, hand-derived expectations).
- **Fixed two correctness bugs** in the HKSJ pooling: (1) the HKSJ variance had
  no floor, so when `Q < k−1` the interval narrowed below the random-effects
  interval — now floored at the RE variance; (2) the t critical value was a crude
  hardcoded approximation — now an exact `qt(0.975, k−1)` table.
- Added Pages scaffold (`.nojekyll`, README, `.gitignore`); renamed
  `culpritcotro.html` → `index.html`.

## Body (E156 draft — CURRENT BODY)

Does treating non-culprit coronary lesions, rather than the culprit alone, change
the balance of recurrent infarction against death and renal harm? This dashboard
reconstructs odds ratios from culprit-only versus complete-revascularization
trials spanning stable STEMI, NSTEMI and cardiogenic-shock populations. It pools
log odds ratios with a Hartung–Knapp–Sidik–Jonkman random-effects model,
reporting τ², I² and a t-based interval, then projects a competing-risk net
clinical benefit through a Monte-Carlo panel. Under stable-trial presets the
efficacy signal favours completeness, but adding the shock trial inflates
heterogeneity and erases the safety margin, so the pooled direction depends on
which populations are admitted. A revival audit found and fixed an
anti-conservative HKSJ interval and a crude t-value, then locked the core behind
a hand-derived 28-assertion test suite. The honest read is that completeness
helps in stable disease but the paradox is population selection, not a universal
rule. The tool is a transparent synthesis aid, not a clinical decision rule.

SUBMITTED: [ ]
