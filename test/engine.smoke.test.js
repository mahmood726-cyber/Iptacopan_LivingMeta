/*
 * Minimal dependency-free smoke test for the Iptacopan LivingMeta stat engine.
 *
 * The dashboard ships its meta-analysis engine inline inside IPTACOPAN_REVIEW.html
 * (functions binaryEffect / hrFromPublished / poolWith / pauleMandelTau2). There is
 * no module boundary to import in a browser-only single-file app, so this test
 * re-implements the SAME core formulas verbatim and asserts they reproduce
 * independently hand-computed reference values. If the engine math is ever changed
 * in the HTML, this file is the canonical reference to check it against.
 *
 * Run: node test/engine.smoke.test.js   (exit 0 = pass, exit 1 = fail)
 */
'use strict';

var Z95 = 1.959964;
var failures = 0;
function approx(name, got, want, tol) {
  tol = tol == null ? 1e-6 : tol;
  if (!Number.isFinite(got) || Math.abs(got - want) > tol) {
    console.error('FAIL ' + name + ': got ' + got + ' want ' + want + ' (tol ' + tol + ')');
    failures++;
  } else {
    console.log('ok   ' + name + ' = ' + got);
  }
}
function ok(name, cond) {
  if (cond) { console.log('ok   ' + name); } else { console.error('FAIL ' + name); failures++; }
}

// --- Engine formulas (verbatim from IPTACOPAN_REVIEW.html) ---
function binaryEffect(tE, tN, cE, cN, measure) {
  if (tN <= 0 || cN <= 0) return null;
  var tNonE = tN - tE, cNonE = cN - cE;
  if (tNonE < 0 || cNonE < 0) return null;
  // 0.5 continuity correction ONLY when >=1 cell is zero (advanced-stats rule)
  var add = (tE === 0 || tNonE === 0 || cE === 0 || cNonE === 0) ? 0.5 : 0;
  var a = tE + add, b = tNonE + add, c = cE + add, e = cNonE + add;
  if (a <= 0 || b <= 0 || c <= 0 || e <= 0) return null;
  var y, v;
  if (measure === 'RR') {
    var pt = a / (a + b), pc = c / (c + e);
    if (pt <= 0 || pc <= 0) return null;
    y = Math.log(pt / pc);
    v = (1 / a) - (1 / (a + b)) + (1 / c) - (1 / (c + e));
  } else {
    y = Math.log(a * e / (b * c));        // logOR
    v = 1 / a + 1 / b + 1 / c + 1 / e;
  }
  if (!Number.isFinite(y) || !Number.isFinite(v) || v <= 0) return null;
  return { y: y, v: v };
}

function hrFromPublished(hr, lo, hi) {
  if (hr <= 0 || lo <= 0 || hi <= 0) return null;
  var y = Math.log(hr);
  var se = (Math.log(hi) - Math.log(lo)) / (2 * Z95);
  if (!Number.isFinite(se) || se <= 0) return null;
  return { y: y, v: se * se };
}

function poolWith(y, v, tau2) {
  var w = v.map(function (vi) { return 1 / (vi + tau2); });
  var sumW = w.reduce(function (a, b) { return a + b; }, 0);
  var effect = 0;
  for (var i = 0; i < y.length; i++) effect += y[i] * w[i];
  effect /= sumW;
  var se = 1 / Math.sqrt(sumW);
  return { effect: effect, se: se, lo: effect - Z95 * se, hi: effect + Z95 * se };
}

function pauleMandelTau2(y, v) {
  var k = y.length;
  if (k < 2) return 0;          // k<1 / k<2 guard: undefined heterogeneity
  var tau2 = 0, target = k - 1;
  for (var iter = 0; iter < 200; iter++) {
    var w = v.map(function (vi) { return 1 / (vi + tau2); });
    var sumW = w.reduce(function (a, b) { return a + b; }, 0);
    if (!(sumW > 0)) break;
    var yBar = 0, i;
    for (i = 0; i < k; i++) yBar += y[i] * w[i];
    yBar /= sumW;
    var Q = 0;
    for (i = 0; i < k; i++) { var d = y[i] - yBar; Q += w[i] * d * d; }
    var diff = Q - target;
    if (Math.abs(diff) < 1e-5) break;
    var slope = 0;
    for (i = 0; i < k; i++) { var dm = y[i] - yBar; slope += -w[i] * w[i] * dm * dm; }
    if (!Number.isFinite(slope) || Math.abs(slope) < 1e-14) break;
    var t2new = Math.max(0, tau2 + diff / slope);
    if (!Number.isFinite(t2new)) break;
    if (Math.abs(t2new - tau2) < 1e-9) { tau2 = t2new; break; }
    tau2 = t2new;
  }
  return tau2;
}

// --- Tests ---

// 1. Per-study log-OR + variance (hand-checked).
var e1 = binaryEffect(10, 100, 20, 100, 'OR');
approx('logOR study1', e1.y, -0.810930, 1e-5);
approx('var study1', e1.v, 0.173611, 1e-5);

// 2. Fixed-effect (tau2=0) inverse-variance pool of two studies, back-transformed.
var e2 = binaryEffect(5, 80, 12, 80, 'OR');
var p = poolWith([e1.y, e2.y], [e1.v, e2.v], 0);
approx('pooled logOR', p.effect, -0.869108, 1e-5);
approx('pooled OR', Math.exp(p.effect), 0.419326, 1e-5);
approx('pooled OR lo', Math.exp(p.lo), 0.217956, 1e-5);
approx('pooled OR hi', Math.exp(p.hi), 0.806741, 1e-5);

// 3. Zero-cell: 0.5 correction applied -> finite, defined effect.
var z = binaryEffect(0, 50, 8, 50, 'OR');
ok('zero-cell yields finite effect', z && Number.isFinite(z.y) && Number.isFinite(z.v));
// And NO correction when all cells are non-zero (a=10 => log(...) uses raw counts).
ok('no continuity correction when cells non-zero', Math.abs(binaryEffect(10, 100, 20, 100, 'OR').y - Math.log(10 * 80 / (90 * 20))) < 1e-12);

// 4. HR from published CI: log scale, SE from CI width.
var h = hrFromPublished(0.5, 0.30, 0.83);
approx('logHR', h.y, Math.log(0.5), 1e-12);
ok('HR variance positive', h.v > 0);

// 5. tau2 guards: k<2 returns 0 (heterogeneity undefined for a single study).
approx('tau2 k=1 is 0', pauleMandelTau2([0.1], [0.2]), 0, 0);
ok('tau2 non-negative', pauleMandelTau2([0.1, -0.4, 0.2, 0.05], [0.2, 0.3, 0.1, 0.25]) >= 0);

// 6. Degenerate input guards return null rather than NaN.
ok('zero N returns null', binaryEffect(5, 0, 3, 50, 'OR') === null);
ok('events exceed N returns null', binaryEffect(60, 50, 3, 50, 'OR') === null);
ok('bad HR CI returns null', hrFromPublished(0, 0.3, 0.8) === null);

if (failures) {
  console.error('\n' + failures + ' assertion(s) FAILED');
  process.exit(1);
}
console.log('\nAll engine smoke assertions passed.');
