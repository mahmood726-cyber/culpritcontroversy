/*
 * RevascLogic engine — pure meta-analysis core for the culprit-only vs complete
 * revascularization evidence-synthesis dashboard ("The Completeness Paradox").
 *
 * Extracted from the dashboard's inline web-worker so the statistical core is a
 * single source of truth, importable under Node for testing and used by the page
 * synchronously (the worker plumbing was removed — the computation is small).
 * Browser: functions are globals (plain declarations). Node: module.exports.
 *
 * Method: Hartung–Knapp–Sidik–Jonkman (HKSJ) random-effects pooling of log odds
 * ratios. τ² and I² via DerSimonian–Laird method-of-moments on fixed-effect
 * weights; HKSJ variance with the standard floor at the classical RE variance;
 * t_{k-1} critical value.
 *
 * Two correctness fixes applied during the 2026-06 revival (see runHKSJ):
 *   1. HKSJ variance FLOOR — the original computed q_term/((k-1)*swRE) with no
 *      floor, so when Q < k-1 (HK multiplier < 1) the HK interval narrowed
 *      *below* the classical random-effects interval (anti-conservative). It is
 *      now floored: var_HKSJ = max(1, qStat) * var_RE, where qStat is the HK
 *      multiplier and var_RE = 1/ΣwRE.
 *   2. Proper t critical value — the original hardcoded crude approximations
 *      (k<10 -> 2.3, etc.). Now uses qt(0.975, k-1) from an exact table
 *      (e.g. k=2 -> 12.706, k=3 -> 4.303, k=5 -> 2.776).
 */

// Inverse Student-t at the two-sided 0.975 level. Exact table for df 1..30
// (matches R qt() to 3 dp); Cornish–Fisher expansion for df > 30.
function tInv975(df) {
    const T = {
        1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
        8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145,
        15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086, 21: 2.080,
        22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060, 26: 2.056, 27: 2.052, 28: 2.048,
        29: 2.045, 30: 2.042
    };
    if (df <= 0) return Infinity;
    const d = Math.round(df);
    if (d <= 30) return T[d];
    const z = 1.959963985;             // qnorm(0.975)
    const z2 = z * z, z3 = z2 * z, z5 = z3 * z2;
    return z
        + (z3 + z) / (4 * df)
        + (5 * z5 + 16 * z3 + 3 * z) / (96 * df * df);
}

// HKSJ pooling on the supplied effect/SE keys.
// data: [{ [kEff]: logOR, [kSe]: seLog, id, active }]
// Returns { mu, se_pool, tVal, tau2, I2, trials } or null on empty input.
function runHKSJ(data, kEff, kSe) {
    if (!data || data.length === 0) return null;

    const k = data.length;
    const y = data.map(d => d[kEff]);
    const v = data.map(d => d[kSe] ** 2);

    const perTrial = data.map(d => ({
        id: d.id,
        or: Math.exp(d[kEff]),
        lo: Math.exp(d[kEff] - 1.96 * d[kSe]),
        hi: Math.exp(d[kEff] + 1.96 * d[kSe]),
        active: d.active
    }));

    if (k === 1) {
        return {
            mu: y[0],
            se_pool: data[0][kSe],
            tVal: 1.96,
            tau2: 0,
            I2: 0,
            trials: perTrial
        };
    }

    // 1. DerSimonian–Laird tau^2 and I^2 (fixed-effect weights for Q and C)
    const wFE = v.map(x => 1 / x);
    const swFE = wFE.reduce((a, b) => a + b, 0);
    const muFE = wFE.reduce((a, b, j) => a + b * y[j], 0) / swFE;

    let Q = 0;
    y.forEach((val, i) => Q += wFE[i] * (val - muFE) ** 2);
    const df = k - 1;
    const C = swFE - wFE.reduce((a, b) => a + b * b, 0) / swFE;
    const tau2 = Math.max(0, (Q - df) / C);
    const I2 = Math.max(0, (Q - df) / Q) * 100;

    // 2. Random-effects weighted mean
    const wRE = v.map(x => 1 / (x + tau2));
    const swRE = wRE.reduce((a, b) => a + b, 0);
    const muRE = wRE.reduce((a, b, j) => a + b * y[j], 0) / swRE;

    // 3. Hartung–Knapp variance with the standard floor at the RE variance.
    let qTerm = 0;
    y.forEach((val, i) => qTerm += wRE[i] * (val - muRE) ** 2);
    const qStat = qTerm / df;                    // HK multiplier
    const varRE = 1 / swRE;                       // classical RE variance
    const varHKSJ = Math.max(1, qStat) * varRE;   // <-- floor (fix #1)
    const seRE = Math.sqrt(varHKSJ);

    const tVal = tInv975(df);                     // <-- proper t_{k-1} (fix #2)

    return {
        mu: muRE,
        se_pool: seRE,
        tVal: tVal,
        tau2: tau2,
        I2: I2,
        trials: perTrial
    };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { runHKSJ, tInv975 };
}
