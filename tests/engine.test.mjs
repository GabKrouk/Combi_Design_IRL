// CombiDesign — engine tests:  node --test tests/
import test from "node:test";
import assert from "node:assert/strict";
import {
  rng, choose, combos, pairCoveringSize, pairCovering, tripletCovering, wavePlan, missingCombinations,
  fSf, tP2, betai, recoveryRate, fitN, formulaN95, formulaValid, oneByOne,
} from "../web/engine.js";

test("combinatorics", () => {
  assert.equal(choose(18, 2), 153);
  assert.equal(choose(18, 3), 816);
  assert.equal(combos(5, 3).length, 10);
  assert.equal(oneByOne(18, 2), 172);
  assert.equal(oneByOne(18, 3), 988);
});

test("minimal pair design sizes (Kleitman & Spencer)", () => {
  // known optimal sizes of binary strength-2 covering arrays
  const known = { 2: 4, 3: 4, 4: 5, 5: 6, 10: 6, 11: 7, 15: 7, 16: 8, 18: 8, 35: 8, 36: 9, 56: 9, 57: 10, 126: 10 };
  for (const [K, N] of Object.entries(known)) assert.equal(pairCoveringSize(+K), N, `K=${K}`);
});

test("pair design covers every pair, first mixture is the control", () => {
  for (const K of [3, 7, 18, 30]) {
    const A = pairCovering(K);
    assert.equal(A.length, pairCoveringSize(K));
    assert.equal(missingCombinations(A, 2), 0);
    assert.equal(A[0].reduce((a, b) => a + b, 0), 0);
  }
});

test("triplet design for 18 compounds covers every triplet with ≤ 22 mixtures", () => {
  const T = tripletCovering(18, { seed: 1 });
  assert.equal(missingCombinations(T, 3), 0);
  assert.ok(T.length <= 22, `found ${T.length}`);
});

test("wave plan: triplet design, then foldover", () => {
  const T = tripletCovering(10, { seed: 2 });
  const { A, wave } = wavePlan(10, 2 * T.length + 5, { base: T });
  assert.equal(A.length, 2 * T.length + 5);
  assert.equal(missingCombinations(A, 3), 0);
  for (let i = 0; i < T.length; i++) assert.deepEqual(A[T.length + i], A[i].map(x => 1 - x));
  assert.deepEqual([...new Set(wave)], [1, 2, 3]);
});

test("statistics match scipy", () => {
  // scipy.stats.f.sf(4, 1, 30); 2 * scipy.stats.t.sf(2.042, 30); scipy.special.betainc(2, 3, 0.4)
  assert.ok(Math.abs(fSf(4, 1, 30) - 0.054625044962983094) < 1e-10);
  assert.ok(Math.abs(tP2(2.042, 30) - 0.050028670656197906) < 1e-10);
  assert.ok(Math.abs(betai(2, 3, 0.4) - 0.5248) < 1e-12);
});

test("seeded random numbers are reproducible", () => {
  const a = rng(5), b = rng(5);
  for (let i = 0; i < 5; i++) assert.equal(a(), b());
  const r = rng(9); let s = 0, s2 = 0; for (let i = 0; i < 20000; i++) { const z = r.normal(); s += z; s2 += z * z; }
  assert.ok(Math.abs(s / 20000) < 0.03 && Math.abs(s2 / 20000 - 1) < 0.05);
});

test("formulas", () => {
  assert.ok(Math.abs(formulaN95(14, 5, 2) - 73.56) < 0.01);
  assert.ok(Math.abs(formulaN95(14, 5, 3) - 103.48) < 0.01);
  assert.equal(formulaValid({ K: 18, reps: 3 }), true);
  assert.equal(formulaValid({ K: 20 }), false);
});

test("simulation agrees with the formula (pairs, S = 14, M = 5)", () => {
  // formula: 74 mixtures for 95 %. With 150 simulations the rate must be well above 80 % at 74 and low at 40.
  const hi = recoveryRate(74, { S: 14, M: 5, order: 2 }, 150, 11);
  const lo = recoveryRate(40, { S: 14, M: 5, order: 2 }, 150, 12);
  assert.ok(hi.p > 0.85, `P(74) = ${hi.p}`);
  assert.ok(lo.p < 0.3, `P(40) = ${lo.p}`);
});

test("logistic fit returns the 95 % point", () => {
  const pts = [30, 40, 50, 60, 80].map(N => ({ N, nsim: 1e4, p: 1 / (1 + Math.exp(-(-40 + 10 * Math.log(N)))) }));
  const f = fitN(pts, 0.95);
  assert.ok(Math.abs(f.N - Math.exp((Math.log(19) + 40) / 10)) < 0.5);
});
