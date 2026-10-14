// CombiDesign — user interface
import { formulaN95, formulaValid, oneByOne, choose, pairCoveringSize, fitN, DEFAULTS } from "./engine.js";

const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmt = x => Math.round(x).toLocaleString("en-US").replace(/,/g, " ");
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const PAL2 = ["#6C7A89", "#E0A030", "#2F6DB5", "#D1495B"];
const PAL3 = ["#6C7A89", "#E0A030", "#2F6DB5", "#7B4FA0", "#1F9E89", "#C4652A", "#D16BA5", "#D1495B"];
const tint = (h, a = 0.3) => "#" + [1, 3, 5].map(i => Math.round(255 - (255 - parseInt(h.substr(i, 2), 16)) * a).toString(16).padStart(2, "0")).join("");

// ------------------------------------------------------------------ theme + banner
const root = document.documentElement;
try { const t = localStorage.getItem("combidesign-theme"); if (t) root.dataset.theme = t; } catch {}
$("#theme").addEventListener("click", () => {
  const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
  root.dataset.theme = dark ? "light" : "dark";
  try { localStorage.setItem("combidesign-theme", root.dataset.theme); } catch {}
  if (SIM) drawSim(); if (DES) drawDesign();
});
(function banner() {             // a mosaic of mixtures: compounds present (bright) or absent
  const svg = $("#mosaic"), w = 1400, h = 140, c = 14, s = [];
  svg.setAttribute("viewBox", `0 0 ${w} ${h}`); svg.setAttribute("preserveAspectRatio", "xMidYMid slice");
  for (let x = 0; x < w / c; x++) for (let y = 0; y < h / c; y++) {
    const on = Math.random() < 0.5, fade = 0.15 + 0.6 * (x / (w / c));
    s.push(`<rect x="${x * c + 1}" y="${y * c + 1}" width="${c - 2}" height="${c - 2}" rx="1.5" fill="${on ? "#46b393" : "#9fb0aa"}" opacity="${(on ? 0.55 : 0.07) * fade}"/>`);
  }
  svg.innerHTML = s.join("");
})();

// ------------------------------------------------------------------ calculator
const num = id => { const v = parseFloat($("#" + id).value); return Number.isFinite(v) ? v : NaN; };
function params() {
  return {
    K: Math.round(num("K")), order: +$("input[name=order]:checked").value, S: Math.round(num("S")), M: Math.round(num("M")),
    reps: Math.round(num("reps")), sd: num("sd"), mainLo: num("mainLo"), mainHi: num("mainHi"), intLo: num("intLo"), intHi: num("intHi"),
  };
}
function check(p) {
  if (!(p.K >= 3 && p.K <= 40)) return "Between 3 and 40 compounds.";
  if (!(p.S >= 0 && p.S <= p.K)) return "The number of compounds with an effect must be between 0 and the number of compounds.";
  const nInt = choose(p.K, p.order);
  if (!(p.M >= 0 && p.M <= Math.min(nInt, 30))) return `Between 0 and ${Math.min(nInt, 30)} combined effects.`;
  if (!(p.reps >= 1)) return "At least 1 replicate.";
  if (!(p.sd > 0) || !(p.mainLo >= 0 && p.mainHi >= p.mainLo) || !(p.intLo >= 0 && p.intHi >= p.intLo)) return "Check the effect sizes and noise.";
  return null;
}
function card(lab, big, sub, cls = "") { return `<div class="card ${cls}"><div class="lab">${lab}</div><div class="big">${big}</div><div class="sub">${sub}</div></div>`; }
function updateCards() {
  const p = params(), e = check(p);
  $("#Mlab").textContent = p.order === 3 ? "combined effects (triplets)" : "combined effects (pairs)";
  if (e) { $("#cards").innerHTML = `<p class="err">${esc(e)}</p>`; $("#flag").textContent = ""; return; }
  const one = oneByOne(p.K, p.order), minDes = p.order === 2 ? pairCoveringSize(p.K) : null;
  let h = "";
  if (formulaValid(p)) {
    const n = formulaN95(p.S, p.M, p.order), rec = Math.ceil(n * 1.1);
    h += card("Mixtures needed", `${rec}<small>mixtures</small>`, `${fmt(rec * p.reps)} tubes with ${p.reps} replicates · formula: ${n.toFixed(0)} for 95 %, + 10 % margin`);
    h += card("One combination at a time", `${fmt(one)}<small>mixtures</small>`, `control + ${p.K} single compounds + ${fmt(choose(p.K, 2))} pairs${p.order === 3 ? ` + ${fmt(choose(p.K, 3))} triplets` : ""}`, "dim");
    h += card("Saving", `÷ ${(one / rec).toFixed(1)}`, `${fmt(one - rec)} mixtures fewer (${fmt((one - rec) * p.reps)} tubes)`);
    $("#flag").textContent = "";
  } else {
    h += card("Mixtures needed", `?`, `the formula is calibrated for 18 compounds, 3 replicates and the default effect sizes: simulate your situation below`, "dim");
    h += card("One combination at a time", `${fmt(one)}<small>mixtures</small>`, `control + ${p.K} single compounds + ${fmt(choose(p.K, 2))} pairs${p.order === 3 ? ` + ${fmt(choose(p.K, 3))} triplets` : ""}`, "dim");
    if (minDes) h += card("Minimal pair design", `${minDes}<small>mixtures</small>`, `shows every pair in its 4 situations, but only enough when few compounds act`);
    $("#flag").textContent = "";
  }
  $("#cards").innerHTML = h;
}
$$("#cform input, #cform select").forEach(el => el.addEventListener("input", updateCards));
updateCards();

// ------------------------------------------------------------------ worker
let worker, JOB = 0;
function getWorker() { if (!worker) worker = new Worker(new URL("worker.js", import.meta.url), { type: "module" }); return worker; }

// ------------------------------------------------------------------ simulation
let SIM = null;
$("#simulate").addEventListener("click", () => {
  const p = params(), e = check(p);
  if (e) { $("#simstatus").innerHTML = `<span class="err">${esc(e)}</span>`; return; }
  const opts = { nsim: +$("#nsim").value, seed: Math.max(1, Math.round(num("seed")) || 1) };
  const id = ++JOB, w = getWorker(), pts = [];
  SIM = { p, opts, points: pts, N95: NaN };
  $("#simulate").disabled = true; $("#simres").hidden = false; $("#simprog").hidden = false;
  $("#simstatus").innerHTML = `<span class="loader"></span>simulating…`;
  $("#simsummary").textContent = ""; drawSim();
  w.onmessage = ({ data: m }) => {
    if (m.id !== id) return;
    if (m.type === "point") {
      pts.push(m.point);
      $("#simprog div").style.width = Math.min(100, 100 * Math.max(...pts.map(q => q.p)) / 0.95) + "%";
      $("#simstatus").innerHTML = `<span class="loader"></span>${m.point.N} mixtures: ${(100 * m.point.p).toFixed(0)} % of the simulated experiments found every effect…`;
      drawSim();
    } else if (m.type === "done") {
      SIM.N95 = m.result.N95; SIM.fit = m.result.fit;
      $("#simulate").disabled = false; $("#simprog").hidden = true; $("#simstatus").textContent = "";
      const rec = Math.ceil(SIM.N95 * 1.1);
      const one = oneByOne(p.K, p.order);
      $("#simsummary").innerHTML = Number.isFinite(SIM.N95)
        ? `With your settings, <b>${Math.ceil(SIM.N95)} mixtures</b> find every effect in 95 % of the simulated experiments.
           With a 10 % safety margin: <b>${rec} mixtures</b>, ${fmt(rec * p.reps)} tubes
           (one combination at a time: ${fmt(one)} mixtures, ÷ ${(one / rec).toFixed(1)}).`
        : `95 % was not reached up to ${pts.at(-1).N} mixtures: the effects are too small for the noise. Add replicates or reduce the noise.`;
      if (Number.isFinite(SIM.N95)) $("#DN").value = rec;
      drawSim();
    } else if (m.type === "error") {
      $("#simulate").disabled = false; $("#simprog").hidden = true; $("#simstatus").innerHTML = `<span class="err">${esc(m.message)}</span>`;
    }
  };
  w.postMessage({ id, cmd: "n95", params: p, opts });
});

function simSvg() {
  const pts = SIM.points, W = 860, H = 330, L = 56, R = 18, T = 14, B = 44;
  const maxN = Math.max(60, ...pts.map(q => q.N), Number.isFinite(SIM.N95) ? SIM.N95 * 1.15 : 0), minN = Math.min(...pts.map(q => q.N), maxN) * 0.9;
  const x = n => L + (W - L - R) * (n - minN) / (maxN - minN), y = v => T + (H - T - B) * (1 - v);
  const fg = css("--fg"), sub = css("--sub"), grid = css("--grid"), line = css("--line"), cut = css("--cut");
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" font-family="Source Sans 3, Arial, sans-serif" font-size="12">`;
  for (let v = 0; v <= 1.0001; v += 0.25) s += `<line x1="${L}" x2="${W - R}" y1="${y(v)}" y2="${y(v)}" stroke="${grid}"/><text x="${L - 8}" y="${y(v) + 4}" text-anchor="end" fill="${sub}">${Math.round(v * 100)} %</text>`;
  const step = maxN - minN > 300 ? 100 : maxN - minN > 120 ? 50 : 20;
  for (let n = Math.ceil(minN / step) * step; n <= maxN; n += step) s += `<line x1="${x(n)}" x2="${x(n)}" y1="${y(0)}" y2="${y(0) + 4}" stroke="${sub}"/><text x="${x(n)}" y="${y(0) + 17}" text-anchor="middle" fill="${sub}">${n}</text>`;
  s += `<text x="${(L + W - R) / 2}" y="${H - 6}" text-anchor="middle" fill="${sub}">number of mixtures</text>`;
  s += `<line x1="${L}" x2="${W - R}" y1="${y(0.95)}" y2="${y(0.95)}" stroke="${cut}" stroke-dasharray="5 4"/><text x="${W - R}" y="${y(0.95) - 5}" text-anchor="end" fill="${cut}">95 %</text>`;
  if (pts.length >= 2) {
    const f = fitN(pts);
    if (f.b > 0) {
      let d = ""; for (let i = 0; i <= 80; i++) { const n = minN + (maxN - minN) * i / 80, v = 1 / (1 + Math.exp(-(f.a + f.b * Math.log(n)))); d += (i ? "L" : "M") + x(n).toFixed(1) + "," + y(v).toFixed(1); }
      s += `<path d="${d}" fill="none" stroke="${line}" stroke-width="2" opacity=".55"/>`;
    }
  }
  for (const q of pts) s += `<circle cx="${x(q.N)}" cy="${y(q.p)}" r="4.5" fill="${line}"><title>${q.N} mixtures: ${(100 * q.p).toFixed(1)} %</title></circle>`;
  if (Number.isFinite(SIM.N95)) s += `<line x1="${x(SIM.N95)}" x2="${x(SIM.N95)}" y1="${y(0)}" y2="${y(0.95)}" stroke="${cut}"/><text x="${x(SIM.N95) + 6}" y="${y(0.08)}" fill="${fg}" font-weight="600">${Math.ceil(SIM.N95)} mixtures</text>`;
  return s + "</svg>";
}
function drawSim() {
  if (!SIM) return;
  $("#simchart").innerHTML = SIM.points.length ? simSvg() : `<p class="note" style="margin:8px">first points coming…</p>`;
  $("#simtab tbody").innerHTML = SIM.points.map(q => `<tr><td class="num">${q.N}</td><td class="num">${q.N * SIM.p.reps}</td><td class="num">${(100 * q.p).toFixed(1)} %</td><td class="num">${(100 * q.mains).toFixed(1)} %</td><td class="num">${(100 * q.ints).toFixed(1)} %</td><td class="num">${q.nsim}</td></tr>`).join("");
}
function download(name, text, type = "text/plain") {
  const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
$("#dlsim").addEventListener("click", () => {
  const p = SIM.p, head = `# CombiDesign simulation: K=${p.K} order=${p.order} S=${p.S} M=${p.M} reps=${p.reps} sd=${p.sd} simple=${p.mainLo}-${p.mainHi} combined=${p.intLo}-${p.intHi} N95=${SIM.N95?.toFixed(1)}\n`;
  download("combidesign_simulation.tsv", head + "mixtures\ttubes\tP_all\tP_simple\tP_combined\tsimulations\n" + SIM.points.map(q => [q.N, q.N * p.reps, q.p, q.mains, q.ints, q.nsim].join("\t")).join("\n") + "\n");
});
$("#dlsimsvg").addEventListener("click", () => download("combidesign_simulation.svg", simSvg(), "image/svg+xml"));
$("#todesign").addEventListener("click", () => { $("input[name=kind][value=waves]").checked = true; $("#design").scrollIntoView({ behavior: "smooth" }); });

// ------------------------------------------------------------------ design
let DES = null;
const names = K => { const l = $("#names").value.split(/\n/).map(s => s.trim()).filter(Boolean); return Array.from({ length: K }, (_, j) => l[j] || `C${j + 1}`); };
$("#dform").addEventListener("submit", ev => {
  ev.preventDefault();
  const K = Math.round(num("K")), kind = $("input[name=kind]:checked").value, N = Math.round(num("DN"));
  if (!(K >= 3 && K <= 40)) { $("#dstatus").innerHTML = `<span class="err">Set between 3 and 40 compounds in the calculator.</span>`; return; }
  if (kind === "waves" && !(N >= 4 && N <= 2000)) { $("#dstatus").innerHTML = `<span class="err">Between 4 and 2000 mixtures.</span>`; return; }
  const id = ++JOB, w = getWorker();
  $("#gen").disabled = true; $("#dstatus").innerHTML = `<span class="loader"></span>building the design…`;
  w.onmessage = ({ data: m }) => {
    if (m.id !== id) return;
    if (m.type === "progress") $("#dstatus").innerHTML = `<span class="loader"></span>searching the triplet design: ${m.N} mixtures…`;
    else if (m.type === "done") {
      $("#gen").disabled = false; $("#dstatus").textContent = "";
      DES = { ...m, K, kind, names: names(K) }; drawDesign(); $("#dres").hidden = false;
    } else if (m.type === "error") { $("#gen").disabled = false; $("#dstatus").innerHTML = `<span class="err">${esc(m.message)}</span>`; }
  };
  w.postMessage({ id, cmd: "design", kind, K, N, seed: Math.max(1, Math.round(num("seed")) || 1) });
});
function hlCols() {
  const v = $("#hl").value.split(/[^0-9]+/).filter(Boolean).map(x => +x - 1);
  const u = [...new Set(v)].filter(j => j >= 0 && j < (DES?.K || 0));
  return u.length === 2 || u.length === 3 ? u : [];
}
$("#hl").addEventListener("input", () => DES && drawDesign());
function heatSvg() {
  const { A, wave, K } = DES, nm = DES.names, N = A.length, H = hlCols(), pal = H.length === 3 ? PAL3 : PAL2;
  const cw = Math.max(14, Math.min(34, Math.floor(900 / K))), ch = N > 120 ? 6 : N > 60 ? 9 : N > 30 ? 13 : 20;
  const L = 52, T = 70, W = L + cw * K + 10, Ht = T + ch * N + 10;
  const on = css("--on"), off = css("--off"), sub = css("--sub"), fg = css("--fg"), cut = css("--cut"), bg = css("--card");
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${Ht}" width="${W}" font-family="Source Sans 3, Arial, sans-serif" font-size="11"><rect width="${W}" height="${Ht}" fill="${bg}"/>`;
  nm.forEach((n, j) => { const x = L + j * cw + cw / 2; s += `<text transform="translate(${x + 4},${T - 6}) rotate(-60)" fill="${H.includes(j) ? fg : sub}" font-weight="${H.includes(j) ? 700 : 400}">${esc(n.length > 14 ? n.slice(0, 13) + "…" : n)}</text>`; });
  const every = N > 60 ? 10 : N > 30 ? 5 : 1;
  for (let i = 0; i < N; i++) {
    if (i % every === 0 || i === N - 1) s += `<text x="${L - 6}" y="${T + i * ch + ch / 2 + 4}" text-anchor="end" fill="${sub}">${i + 1}</text>`;
    let code = 0; for (const j of H) code = code * 2 + A[i][j];
    for (let j = 0; j < K; j++) {
      let f = A[i][j] ? on : off, op = 1;
      if (H.length) { if (H.includes(j)) f = A[i][j] ? pal[code] : tint(pal[code]); else op = 0.35; }
      s += `<rect x="${L + j * cw + 0.5}" y="${T + i * ch + 0.5}" width="${cw - 1}" height="${ch - 1}" fill="${f}" opacity="${op}"/>`;
    }
  }
  if (wave) for (let i = 1; i < N; i++) if (wave[i] !== wave[i - 1]) s += `<line x1="${L - 2}" x2="${L + cw * K + 2}" y1="${T + i * ch}" y2="${T + i * ch}" stroke="${cut}" stroke-width="2"/>`;
  for (const j of H) s += `<rect x="${L + j * cw}" y="${T}" width="${cw}" height="${ch * N}" fill="none" stroke="${fg}" stroke-width="1.2"/>`;
  return s + "</svg>";
}
function drawDesign() {
  const { A, wave, K, check, kind } = DES, N = A.length, per = A.map(r => r.reduce((a, b) => a + b, 0));
  const pairsOk = check.pairs === 0, tripOk = check.triplets === 0;
  let msg = `<b>${N} mixtures</b> for ${K} compounds, ${Math.min(...per)} to ${Math.max(...per)} compounds per mixture. `;
  msg += `Every pair seen in its 4 situations: <b>${pairsOk ? "yes" : "no (" + check.pairs + " missing)"}</b>`;
  if (check.triplets !== null) msg += ` · every triplet in its 8 situations: <b>${tripOk ? "yes" : "no (" + check.triplets + " missing)"}</b>`;
  if (wave) { const nw = Math.max(...wave); msg += `. ${nw} wave${nw > 1 ? "s" : ""} of ${DES.A.length >= 1 ? wave.filter(v => v === 1).length : 0} mixtures (red lines); run them in order.`; }
  else msg += ".";
  $("#dsummary").innerHTML = msg;
  const H = hlCols();
  $("#dlegend").innerHTML = `<span><i style="background:var(--on)"></i>compound present</span><span><i style="background:var(--off);border:1px solid var(--rule)"></i>compound absent</span>${wave ? `<span><i style="background:var(--cut);height:3px;vertical-align:3px"></i>new wave</span>` : ""}`;
  if (H.length) {
    const pal = H.length === 3 ? PAL3 : PAL2, k = H.length;
    let h = "";
    for (let c = 0; c < (1 << k); c++) {
      const bits = H.map((_, b) => (c >> (k - 1 - b)) & 1), rows = [];
      A.forEach((r, i) => { let v = 0; for (const j of H) v = v * 2 + r[j]; if (v === c) rows.push(i + 1); });
      const pres = H.filter((_, b) => bits[b]).map(j => DES.names[j]);
      h += `<div class="chips">${H.map((j, b) => `<span style="background:${bits[b] ? pal[c] : tint(pal[c])};color:${bits[b] ? "#fff" : "#444"}">${j + 1}</span>`).join("")}</div>
        <div><b>${pres.length ? esc(pres.join(" + ")) : "none"}</b> → ${rows.length ? `mixtures ${rows.slice(0, 12).join(", ")}${rows.length > 12 ? "…" : ""} (${rows.length})` : `<span class="err">never</span>`}</div>`;
    }
    $("#hllegend").innerHTML = h;
  } else $("#hllegend").innerHTML = "";
  $("#heat").innerHTML = heatSvg();
}
$("#dlcsv").addEventListener("click", () => {
  const { A, wave } = DES, nm = DES.names;
  const head = ["mixture", ...(wave ? ["wave"] : []), ...nm, "n_compounds"].map(s => /[",;]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s).join(",");
  download("combidesign_design.csv", head + "\n" + A.map((r, i) => [i + 1, ...(wave ? [wave[i]] : []), ...r, r.reduce((a, b) => a + b, 0)].join(",")).join("\n") + "\n", "text/csv");
});
$("#dlheat").addEventListener("click", () => download("combidesign_design.svg", heatSvg(), "image/svg+xml"));
