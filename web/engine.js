// CombiDesign — engine: combinatorial designs, Monte Carlo simulation of effect recovery, statistics.
// Pure ES module, no dependency: runs in the browser (main thread or worker) and in Node (tests).

// ------------------------------------------------------------------ random numbers (seeded)
export function rng(seed = 1) {
  let a = seed >>> 0 || 1;
  const u = () => {                       // mulberry32
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  let spare = null;
  u.normal = () => {                      // Box–Muller
    if (spare !== null) { const s = spare; spare = null; return s; }
    let x, y, r;
    do { x = 2 * u() - 1; y = 2 * u() - 1; r = x * x + y * y; } while (r === 0 || r >= 1);
    const f = Math.sqrt(-2 * Math.log(r) / r);
    spare = y * f; return x * f;
  };
  u.int = n => Math.floor(u() * n);
  u.shuffle = arr => { for (let i = arr.length - 1; i > 0; i--) { const j = u.int(i + 1); [arr[i], arr[j]] = [arr[j], arr[i]]; } return arr; };
  return u;
}

// ------------------------------------------------------------------ combinatorics
export function choose(n, k) {
  if (k < 0 || k > n) return 0;
  k = Math.min(k, n - k);
  let r = 1;
  for (let i = 1; i <= k; i++) r = r * (n - k + i) / i;
  return Math.round(r);
}
export function combos(n, k) {
  const out = [], c = Array.from({ length: k }, (_, i) => i);
  if (k > n) return out;
  for (;;) {
    out.push(c.slice());
    let i = k - 1;
    while (i >= 0 && c[i] === n - k + i) i--;
    if (i < 0) return out;
    c[i]++;
    for (let j = i + 1; j < k; j++) c[j] = c[j - 1] + 1;
  }
}

// ------------------------------------------------------------------ designs
/** Number of factors (VOCs) per row etc. A design is an array of rows, each row an array of 0/1 (length K). */

/** Missing (tuple, level combination) pairs for strength t. 0 means every t-subset of factors shows all 2^t combinations. */
export function missingCombinations(A, t) {
  const K = A[0].length;
  let miss = 0;
  for (const c of combos(K, t)) {
    const seen = new Uint8Array(1 << t);
    for (const row of A) { let v = 0; for (const j of c) v = v * 2 + row[j]; seen[v] = 1; }
    for (const s of seen) if (!s) miss++;
  }
  return miss;
}

/** Smallest N of a binary covering array of strength 2 for K factors (Kleitman & Spencer 1973; Katona 1973). */
export function pairCoveringSize(K) {
  if (K < 2) return 1;
  let N = 2;
  while (choose(N - 1, Math.ceil(N / 2)) < K) N++;
  return N;
}

/** Optimal strength-2 design: row 1 = control (no VOC); each column is a distinct ⌈N/2⌉-subset of the other N−1 rows.
 *  Two such subsets always intersect (11) and are never nested (10 and 01): every pair is covered. */
export function pairCovering(K, seed = 1) {
  const N = pairCoveringSize(K), w = Math.ceil(N / 2), R = rng(seed);
  const all = combos(N - 1, w);
  let best = null, bestSpread = Infinity;
  for (let tr = 0; tr < 300; tr++) {           // pick K subsets that balance the number of VOCs per mixture
    const pick = R.shuffle(all.slice()).slice(0, K), load = new Array(N - 1).fill(0);
    for (const s of pick) for (const r of s) load[r]++;
    const spread = Math.max(...load) - Math.min(...load);
    if (spread < bestSpread) { bestSpread = spread; best = pick; if (spread <= 1) break; }
  }
  best.sort((a, b) => a.join() < b.join() ? -1 : 1);
  const A = Array.from({ length: N }, () => new Array(K).fill(0));
  best.forEach((s, j) => { for (const r of s) A[r + 1][j] = 1; });
  return A;
}

/** Strength-3 binary covering design found by greedy construction + targeted simulated annealing
 *  (an AETG-like heuristic, as in Shasha's CombGeneric). Returns the smallest design found. */
export function tripletCovering(K, { seed = 1, iters = 40000, onProgress } = {}) {
  const R = rng(seed);
  const T = combos(K, 3), nT = T.length;
  const byCol = Array.from({ length: K }, () => []);
  T.forEach((t, i) => t.forEach(j => byCol[j].push(i)));
  const code = (row, t) => (row[t[0]] << 2) | (row[t[1]] << 1) | row[t[2]];

  // greedy start
  let A = [];
  {
    const cov = new Uint8Array(nT * 8); let left = nT * 8;
    while (left > 0) {
      let best = null, bestGain = -1;
      for (let c = 0; c < 25; c++) {
        const row = Array.from({ length: K }, () => R.int(2));
        let start = R.int(nT * 8);
        for (let s = 0; s < nT * 8; s++) { const x = (start + s) % (nT * 8); if (!cov[x]) { start = x; break; } }
        const t = T[Math.floor(start / 8)], v = start % 8;
        row[t[0]] = (v >> 2) & 1; row[t[1]] = (v >> 1) & 1; row[t[2]] = v & 1;
        let g = 0; for (let i = 0; i < nT; i++) if (!cov[i * 8 + code(row, T[i])]) g++;
        if (g > bestGain) { bestGain = g; best = row; }
      }
      A.push(best);
      for (let i = 0; i < nT; i++) { const x = i * 8 + code(best, T[i]); if (!cov[x]) { cov[x] = 1; left--; } }
    }
  }
  if (onProgress) onProgress({ N: A.length, phase: "greedy" });

  // annealing: try to repair a design with one row less, repeatedly
  const anneal = A0 => {
    const N = A0.length, M = A0.map(r => r.slice()), cnt = new Int32Array(nT * 8);
    for (const row of M) for (let i = 0; i < nT; i++) cnt[i * 8 + code(row, T[i])]++;
    let miss = 0; for (const c of cnt) if (!c) miss++;
    const stamp = new Int32Array(nT); let st = 0;
    const apply = (r, cols, bits) => {          // returns delta of missing count
      st++; const aff = [];
      for (const j of cols) for (const i of byCol[j]) if (stamp[i] !== st) { stamp[i] = st; aff.push(i); }
      let d = 0;
      for (const i of aff) { const x = i * 8 + code(M[r], T[i]); if (--cnt[x] === 0) d++; }
      const old = cols.map(j => M[r][j]);
      cols.forEach((j, q) => { M[r][j] = bits[q]; });
      for (const i of aff) { const x = i * 8 + code(M[r], T[i]); if (cnt[x]++ === 0) d--; }
      return [d, old];
    };
    let temp = 1;
    for (let it = 0; it < iters && miss > 0; it++) {
      let x = R.int(nT * 8);
      for (let s = 0; s < nT * 8; s++) { const y = (x + s) % (nT * 8); if (!cnt[y]) { x = y; break; } }
      const t = T[Math.floor(x / 8)], v = x % 8, r = R.int(N);
      const [d, old] = apply(r, t, [(v >> 2) & 1, (v >> 1) & 1, v & 1]);
      if (d > 0 && R() >= Math.exp(-d / temp)) apply(r, t, old);
      else miss += d;
      temp = Math.max(0.05, temp * 0.9997);
    }
    return miss === 0 ? M : null;
  };
  for (;;) {
    const drop = R.int(A.length);
    const B = anneal(A.filter((_, i) => i !== drop));
    if (!B) break;
    A = B;
    if (onProgress) onProgress({ N: A.length, phase: "annealing" });
  }
  A.sort((a, b) => a.reduce((s, x) => s + x, 0) - b.reduce((s, x) => s + x, 0));
  return A;
}

/** Design of N mixtures built in waves: strength-3 design, its foldover (presence/absence swapped),
 *  then the same with the VOC columns permuted, and so on; truncated to N rows. Returns {A, wave}. */
export function wavePlan(K, N, { seed = 1, base } = {}) {
  const R = rng(seed + 101);
  const T = base || tripletCovering(K, { seed });
  const A = [], wave = [];
  let w = 0, perm = Array.from({ length: K }, (_, i) => i);
  while (A.length < N) {
    for (const flip of [0, 1]) {
      w++;
      for (const row of T) {
        if (A.length >= N) break;
        A.push(perm.map(j => flip ? 1 - row[j] : row[j])); wave.push(w);
      }
    }
    perm = R.shuffle(perm.slice());
  }
  return { A, wave, waveSize: T.length };
}

/** Random balanced design (each VOC present in ⌊N/2⌋ mixtures): the designs used to calibrate the formulas. */
export function balancedRandom(N, K, R) {
  const A = Array.from({ length: N }, () => new Array(K).fill(0));
  const idx = Array.from({ length: N }, (_, i) => i);
  for (let j = 0; j < K; j++) { R.shuffle(idx); for (let i = 0; i < Math.floor(N / 2); i++) A[idx[i]][j] = 1; }
  return A;
}

// ------------------------------------------------------------------ statistics
function lgamma(x) {                     // Lanczos approximation
  const g = [676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
    12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
  x -= 1; let a = 0.99999999999980993; const t = x + 7.5;
  for (let i = 0; i < 8; i++) a += g[i] / (x + i + 1);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}
function betacf(a, b, x) {
  const MAXIT = 300, EPS = 1e-14, FPMIN = 1e-300;
  let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN; d = 1 / d;
    const del = d * c; h *= del; if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}
/** Regularised incomplete beta I_x(a, b). */
export function betai(a, b, x) {
  if (x <= 0) return 0; if (x >= 1) return 1;
  const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
}
/** Upper tail of the F distribution. */
export const fSf = (F, d1, d2) => F <= 0 ? 1 : betai(d2 / 2, d1 / 2, d2 / (d2 + d1 * F));
/** Two-sided P of Student's t. */
export const tP2 = (t, df) => betai(df / 2, 0.5, df / (df + t * t));

// ------------------------------------------------------------------ simulation
export const DEFAULTS = {
  K: 18, order: 2, S: 14, M: 5, reps: 3, sd: 1,
  mainLo: 0.8, mainHi: 3, intLo: 2.5, intHi: 5,
  alphaMain: 0.05, alphaInt: 0.01,
};

const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };

/** One simulated experiment: N mixtures × reps replicates, random balanced design, S true simple effects,
 *  M true interactions of the given order (2 or 3). Analysis: the K simple effects are always estimated;
 *  interactions are selected forward among all pairs (order 2) or all pairs and triplets (order 3), each step
 *  tested by F with a Bonferroni correction (family-wise risk alphaInt). Success = every true simple effect
 *  significant (t test, alphaMain), and exactly the true interactions selected. */
export function simulateOnce(N, p, R) {
  const { K, order, S, M, reps, sd } = p;
  const A = balancedRandom(N, K, R), n = N * reps;
  const X = Array.from({ length: K }, () => new Float64Array(n));
  for (let i = 0; i < N; i++) for (let r = 0; r < reps; r++) for (let j = 0; j < K; j++) X[j][i * reps + r] = A[i][j];
  const sign = () => R() < 0.5 ? -1 : 1;
  const mains = R.shuffle(Array.from({ length: K }, (_, j) => j)).slice(0, S);
  const tuples = combos(K, order);
  const trueIdx = R.shuffle(Array.from({ length: tuples.length }, (_, i) => i)).slice(0, M);
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) y[i] = 20 + sd * R.normal();
  for (const j of mains) { const b = sign() * (p.mainLo + R() * (p.mainHi - p.mainLo)); for (let i = 0; i < n; i++) y[i] += b * X[j][i]; }
  const trueSet = new Set();
  for (const t of trueIdx) {
    const tu = tuples[t], g = sign() * (p.intLo + R() * (p.intHi - p.intLo)); trueSet.add(tu.join(","));
    for (let i = 0; i < n; i++) { let v = 1; for (const j of tu) v *= X[j][i]; y[i] += g * v; }
  }
  // base: intercept + K simple effects, orthonormalised (modified Gram–Schmidt)
  const Q = [];
  const ortho = v => { const w = Float64Array.from(v); for (const q of Q) { const c = dot(q, w); for (let i = 0; i < n; i++) w[i] -= c * q[i]; } return w; };
  const pushQ = v => { const w = ortho(v), nr = Math.sqrt(dot(w, w)); if (nr < 1e-8) return false; for (let i = 0; i < n; i++) w[i] /= nr; Q.push(w); return true; };
  pushQ(new Float64Array(n).fill(1));
  for (let j = 0; j < K; j++) if (!pushQ(X[j])) return { ok: false, mainsOk: false, intsOk: false };
  // candidates
  const cand = [...combos(K, 2), ...(order >= 3 ? combos(K, 3) : [])];
  const P = cand.length, Z = cand.map(tu => { const v = new Float64Array(n); for (let i = 0; i < n; i++) { let x = 1; for (const j of tu) x *= X[j][i]; v[i] = x; } return ortho(v); });
  let r = ortho(y);
  const alive = new Uint8Array(P).fill(1), sel = [];
  let k = Q.length;
  while (n - k - 2 >= 2 && sel.length <= 3 * M + 10) {
    let best = -1, bs = -1;
    for (let c = 0; c < P; c++) {
      if (!alive[c]) continue;
      const z = Z[c], nz = dot(z, z);
      if (nz < 1e-8 * n) continue;
      const zr = dot(z, r), s = zr * zr / nz;
      if (s > bs) { bs = s; best = c; }
    }
    if (best < 0) break;
    const rss = dot(r, r), df = n - k - 1, F = bs / ((rss - bs) / df);
    if (fSf(F, 1, df) * P >= p.alphaInt) break;
    const z = Z[best], nz = Math.sqrt(dot(z, z)), q = z.map(v => v / nz);
    const qr = dot(q, r); for (let i = 0; i < n; i++) r[i] -= qr * q[i];
    for (let c = 0; c < P; c++) if (alive[c]) { const zc = Z[c], qz = dot(q, zc); if (qz !== 0) for (let i = 0; i < n; i++) zc[i] -= qz * q[i]; }
    alive[best] = 0; sel.push(best); k++;
  }
  const intsOk = sel.length === trueSet.size && sel.every(c => trueSet.has(cand[c].join(",")));
  // final model: t tests of the simple effects
  const cols = [new Float64Array(n).fill(1), ...X, ...sel.map(c => { const v = new Float64Array(n); for (let i = 0; i < n; i++) { let x = 1; for (const j of cand[c]) x *= X[j][i]; v[i] = x; } return v; })];
  const m = cols.length, dfr = n - m;
  if (dfr < 1) return { ok: false, mainsOk: false, intsOk };
  const Qf = [], Rm = Array.from({ length: m }, () => new Float64Array(m));
  for (let a = 0; a < m; a++) {
    const w = Float64Array.from(cols[a]);
    for (let b = 0; b < a; b++) { const c = dot(Qf[b], w); Rm[b][a] = c; for (let i = 0; i < n; i++) w[i] -= c * Qf[b][i]; }
    const nr = Math.sqrt(dot(w, w)); if (nr < 1e-8) return { ok: false, mainsOk: false, intsOk };
    Rm[a][a] = nr; for (let i = 0; i < n; i++) w[i] /= nr; Qf.push(w);
  }
  const qy = Qf.map(q => dot(q, y)), beta = new Float64Array(m);
  for (let a = m - 1; a >= 0; a--) { let s = qy[a]; for (let b = a + 1; b < m; b++) s -= Rm[a][b] * beta[b]; beta[a] = s / Rm[a][a]; }
  const rssF = dot(y, y) - dot(qy, qy), s2 = Math.max(rssF, 1e-12) / dfr;
  const Ri = Array.from({ length: m }, () => new Float64Array(m));      // inverse of the triangular R
  for (let a = m - 1; a >= 0; a--) {
    Ri[a][a] = 1 / Rm[a][a];
    for (let b = a + 1; b < m; b++) { let s = 0; for (let c = a + 1; c <= b; c++) s += Rm[a][c] * Ri[c][b]; Ri[a][b] = -s / Rm[a][a]; }
  }
  let mainsOk = true;
  for (const j of mains) {
    const a = 1 + j; let v = 0; for (let b = a; b < m; b++) v += Ri[a][b] * Ri[a][b];
    if (tP2(beta[a] / Math.sqrt(s2 * v), dfr) >= p.alphaMain) { mainsOk = false; break; }
  }
  return { ok: mainsOk && intsOk, mainsOk, intsOk };
}

/** Probability of full recovery at N mixtures, estimated from nsim simulations. */
export function recoveryRate(N, params, nsim = 200, seed = 1) {
  const p = { ...DEFAULTS, ...params }, R = rng(seed);
  let ok = 0, mo = 0, io = 0;
  for (let s = 0; s < nsim; s++) { const r = simulateOnce(N, p, R); ok += r.ok; mo += r.mainsOk; io += r.intsOk; }
  return { N, p: ok / nsim, mains: mo / nsim, ints: io / nsim, nsim };
}

/** Logistic fit of P(success) against log N (weighted IRLS) and the N where it reaches `target`. */
export function fitN(points, target = 0.95) {
  const xs = points.map(q => Math.log(q.N)), ys = points.map(q => q.p), ws = points.map(q => q.nsim || 200);
  let a = -10, b = 3;
  for (let it = 0; it < 100; it++) {
    let g0 = 0, g1 = 0, h00 = 0, h01 = 0, h11 = 0;
    for (let i = 0; i < xs.length; i++) {
      const pr = 1 / (1 + Math.exp(-(a + b * xs[i]))), w = ws[i] * pr * (1 - pr) + 1e-9, e = ws[i] * (ys[i] - pr);
      g0 += e; g1 += e * xs[i]; h00 += w; h01 += w * xs[i]; h11 += w * xs[i] * xs[i];
    }
    const det = h00 * h11 - h01 * h01; if (Math.abs(det) < 1e-12) break;
    const da = (h11 * g0 - h01 * g1) / det, db = (h00 * g1 - h01 * g0) / det;
    a += da; b += db; if (Math.abs(da) + Math.abs(db) < 1e-9) break;
  }
  const logit = Math.log(target / (1 - target));
  return { a, b, N: b > 0 ? Math.exp((logit - a) / b) : NaN };
}

/** Scan N upward until the recovery rate reaches the target twice, then fit. onPoint(point) after each N. */
export function searchN95(params, { nsim = 200, seed = 1, target = 0.95, maxN = 1000, onPoint } = {}) {
  const p = { ...DEFAULTS, ...params };
  let N = p.K + 1 + p.M + 4, hits = 0;
  const pts = [];
  while (N <= maxN) {
    const pt = recoveryRate(N, p, nsim, seed + N);
    pts.push(pt); if (onPoint) onPoint(pt);
    if (pt.p >= target) hits++;
    if (hits >= 2) break;
    N = Math.max(N + 2, Math.round(N * 1.15));
  }
  const fit = pts.some(q => q.p > 0) && pts.some(q => q.p < 1) ? fitN(pts, target) : { N: pts.length ? pts.find(q => q.p >= target)?.N ?? NaN : NaN };
  return { points: pts, N95: fit.N, fit };
}

// ------------------------------------------------------------------ formulas fitted on the simulations (K = 18)
export const FORMULA = {
  2: { a: 35.06, b: 0.45, c: 6.44, d: 0, r2: 0.996, configs: 54 },
  3: { a: 38.99, b: 0.153, c: 7.94, d: 0.906, r2: 0.996, configs: 42 },
};
/** N95 predicted by the formula (calibrated for 18 VOCs, 3 replicates, noise SD 1, simple effects 0.8–3, interactions 2.5–5). */
export function formulaN95(S, M, order = 2) {
  const f = FORMULA[order];
  return f.a + f.b * S + f.c * M + f.d * M * M;
}
/** True when the parameters are the ones the formula was calibrated on. */
export function formulaValid(p) {
  const q = { ...DEFAULTS, ...p };
  return q.K === 18 && q.reps === 3 && q.sd === 1 && q.mainLo === 0.8 && q.mainHi === 3 && q.intLo === 2.5 && q.intHi === 5;
}
/** Mixtures needed when every combination is tested on its own: control + K single VOCs + all pairs (+ all triplets). */
export function oneByOne(K, order = 2) {
  return 1 + K + choose(K, 2) + (order >= 3 ? choose(K, 3) : 0);
}
