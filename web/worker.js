// CombiDesign — background computations (simulations and design search) so the page stays responsive.
import { searchN95, pairCovering, tripletCovering, wavePlan, missingCombinations } from "./engine.js";

self.onmessage = ({ data: m }) => {
  const post = o => self.postMessage({ id: m.id, ...o });
  try {
    if (m.cmd === "n95") {
      const res = searchN95(m.params, { ...m.opts, onPoint: point => post({ type: "point", point }) });
      post({ type: "done", result: res });
    } else if (m.cmd === "design") {
      const { kind, K, N, seed } = m;
      let A, wave = null, strength;
      if (kind === "pairs") { A = pairCovering(K, seed); strength = 2; }
      else {
        const T = tripletCovering(K, { seed, onProgress: p => post({ type: "progress", ...p }) });
        if (kind === "triplets") { A = T; strength = 3; }
        else { const w = wavePlan(K, N, { seed, base: T }); A = w.A; wave = w.wave; strength = A.length >= T.length ? 3 : null; }
      }
      const check = {
        pairs: A.length ? missingCombinations(A, 2) : null,
        triplets: kind !== "pairs" ? missingCombinations(A, 3) : null,
      };
      post({ type: "done", A, wave, check, strength });
    }
  } catch (e) { post({ type: "error", message: String(e && e.message || e) }); }
};
