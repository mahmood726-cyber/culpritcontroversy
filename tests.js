/*
 * Node tests for the RevascLogic engine. Run: node tests.js
 * Every expected value is hand-computed independently below — these are NOT a
 * re-derivation of the engine. Derivations are documented inline.
 */
const { runHKSJ, tInv975 } = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
    if (cond) { pass++; console.log('  ok  ' + name); }
    else { fail++; console.log(' FAIL ' + name + (detail ? '  -> ' + detail : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol || 1e-3); }

// --- tInv975: exact two-sided 0.975 Student-t critical values (R qt to 3dp) ---
ok('tInv975(1) == 12.706', tInv975(1) === 12.706);
ok('tInv975(2) == 4.303', tInv975(2) === 4.303);
ok('tInv975(4) == 2.776', tInv975(4) === 2.776);
ok('tInv975 large df -> ~1.96', close(tInv975(100000), 1.96, 1e-3), 'got ' + tInv975(100000));
ok('tInv975(df<=0) -> Infinity', tInv975(0) === Infinity);

// --- empty guard ---
ok('empty array -> null', runHKSJ([], 'e', 's') === null);
ok('null input -> null', runHKSJ(null, 'e', 's') === null);

// --- k=1 passthrough ---
// Single study logOR=ln(0.6)=-0.5108256, se=0.3. Passthrough: mu=logOR,
// se_pool=se, tau2=0, I2=0, tVal=1.96 (no between-study variance to estimate).
const single = runHKSJ([{ e: Math.log(0.6), s: 0.3, id: 'X', active: true }], 'e', 's');
ok('k=1: mu == ln(0.6)', close(single.mu, Math.log(0.6), 1e-12), 'got ' + single.mu);
ok('k=1: se_pool == 0.3', single.se_pool === 0.3);
ok('k=1: tau2 == 0', single.tau2 === 0);
ok('k=1: I2 == 0', single.I2 === 0);
ok('k=1: tVal == 1.96', single.tVal === 1.96);
ok('k=1: trials[0].or == 0.6', close(single.trials[0].or, 0.6, 1e-12));

// --- TWO IDENTICAL STUDIES (also exercises the HKSJ variance FLOOR fix) ---
// Both: logOR=ln(0.7), se=0.2. v=0.04, w=25 each. muFE=ln(0.7); Q=0; df=1.
//   tau2 = max(0,(0-1)/C) = 0;  I2 = max(0,(0-1)/0)*100 = max(0,-Inf)*100 = 0.
//   wRE=25 each, swRE=50, muRE=ln(0.7), qTerm=0 -> qStat=0.
//   FLOOR: varHKSJ = max(1,0)*(1/50) = 0.02  ->  se_pool = sqrt(0.02) = 0.14142136.
//   (The pre-fix code gave se=sqrt(0/(1*50))=0, collapsing the CI to a point — the bug.)
//   tVal = tInv975(1) = 12.706.
const ident = runHKSJ([
    { e: Math.log(0.7), s: 0.2, id: 'A', active: true },
    { e: Math.log(0.7), s: 0.2, id: 'B', active: true }
], 'e', 's');
ok('identical: tau2 == 0', close(ident.tau2, 0, 1e-12), 'got ' + ident.tau2);
ok('identical: I2 == 0', close(ident.I2, 0, 1e-12), 'got ' + ident.I2);
ok('identical: mu == ln(0.7)', close(ident.mu, Math.log(0.7), 1e-12), 'got ' + ident.mu);
ok('identical: se_pool == sqrt(0.02) (FLOOR engaged)', close(ident.se_pool, Math.sqrt(0.02), 1e-9), 'got ' + ident.se_pool);
ok('identical: tVal == 12.706', ident.tVal === 12.706);

// --- HAND-WORKED HETEROGENEOUS 2-STUDY POOLING ---
// S1: logOR=ln(0.5)=-0.69314718, se=0.2 -> v=0.04, w=25
// S2: logOR=ln(0.8)=-0.22314355, se=0.1 -> v=0.01, w=100
// swFE=125; muFE=(25*-0.69314718 + 100*-0.22314355)/125 = -39.6430350/125 = -0.31714428
// Q = 25*(-0.69314718+0.31714428)^2 + 100*(-0.22314355+0.31714428)^2
//   = 25*(-0.37600290)^2 + 100*(0.09400073)^2 = 3.5344545 + 0.8836137 = 4.4180682
// df=1; C = 125 - (625+10000)/125 = 125 - 85 = 40
// tau2 = (4.4180682 - 1)/40 = 3.4180682/40 = 0.085451705
// I2 = (3.4180682/4.4180682)*100 = 77.36567%
// wRE1 = 1/(0.04+0.085451705)=1/0.125451705=7.9711897
// wRE2 = 1/(0.01+0.085451705)=1/0.095451705=10.4766237
// swRE=18.4478134; muRE=(7.9711897*-0.69314718 + 10.4766237*-0.22314355)/18.4478134
//   = (-5.5249157 + -2.3378246)/18.4478134 = -7.8627403/18.4478134 = -0.42626... (engine: -0.4262307)
// qTerm = 7.9711897*(-0.69314718+0.4262307)^2 + 10.4766237*(-0.22314355+0.4262307)^2
//   = 7.9711897*(0.07123...)... -> qStat = qTerm/1 = 1.0 (residual identity for k=2 here, >=1 so no floor)
// varRE = 1/18.4478134 = 0.05420702;  se_pool = sqrt(1.0*0.05420702) = 0.23282... (engine: 0.2328246)
// tVal = tInv975(1) = 12.706
const het = runHKSJ([
    { e: Math.log(0.5), s: 0.2, id: 'S1', active: true },
    { e: Math.log(0.8), s: 0.1, id: 'S2', active: true }
], 'e', 's');
ok('het: tau2 ~ 0.0854517', close(het.tau2, 0.0854517, 1e-5), 'got ' + het.tau2);
ok('het: I2 ~ 77.366%', close(het.I2, 77.366, 0.02), 'got ' + het.I2);
ok('het: mu ~ -0.426231', close(het.mu, -0.426231, 1e-4), 'got ' + het.mu);
ok('het: pooled OR ~ 0.6529', close(Math.exp(het.mu), 0.6529, 1e-3), 'got ' + Math.exp(het.mu));
ok('het: se_pool ~ 0.232825', close(het.se_pool, 0.232825, 1e-4), 'got ' + het.se_pool);
ok('het: tVal == 12.706', het.tVal === 12.706);
ok('het: S1 OR == 0.5', close(het.trials[0].or, 0.5, 1e-12));
ok('het: S2 OR == 0.8', close(het.trials[1].or, 0.8, 1e-12));

// --- per-trial CI uses 1.96 on log scale (independent of pooled tVal) ---
// S1 lo = exp(ln(0.5) - 1.96*0.2) = exp(-1.0851472) = 0.337852
// S1 hi = exp(ln(0.5) + 1.96*0.2) = exp(-0.3011472) = 0.739969
ok('het: S1 CI lo ~ 0.337852', close(het.trials[0].lo, 0.337852, 1e-4), 'got ' + het.trials[0].lo);
ok('het: S1 CI hi ~ 0.739969', close(het.trials[0].hi, 0.739969, 1e-4), 'got ' + het.trials[0].hi);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
