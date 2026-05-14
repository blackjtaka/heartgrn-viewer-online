// snp2cell hub viewer.
// Parametric Cytoscape.js viewer driven by manifest.json + per-disease DE-z
// asset + per-(disease × target_cs) payload. UI controls rebuild the subgraph
// and recompute specificity client-side without re-running snp2cell.
//
// State machine (simplified):
//   init()                — load manifest + lineage, wire UI
//   selectDisease(d)      — fetch & cache de_z/{d}.json.gz (~5 MB gz)
//   selectTargetCs(cs)    — fetch & cache payloads/{d}__{cs}.json
//   rebuild()             — recompute subgraph (top_N × peaks × tfs) + spec + layout

(function () {
  "use strict";

  // ---------- module state ----------
  const S = {
    manifest: null,
    lineage: null,
    diseaseId: null,
    targetCs: null,
    scoreKey: null,
    refMode: "all",           // "all" | "selected"
    refCsList: [],            // when refMode === "selected"
    edgeMode: "highlight",    // "highlight" | "subset" | "off"
    edgeHighlightCts: [],
    cellTypePriority: true,   // when true & edgeHighlightCts not empty,
                              // prioritize matching edges in top-N selection
    fdrThresh: 0.10,

    topN: 50,
    peaksPerSeed: 6,
    tfsPerPeak: 4,

    tfLabelMode: "top_specific",
    tfTopN: 10,
    showGeneLabels: true,
    showPeakLabels: false,
    snpEmphasis: true,
    snpHalo: false,

    colorMinOverride: null,
    colorMaxOverride: null,
    colorScaleAuto: { min: 0, max: 1 },
    colorScale: { min: 0, max: 1 },

    // API base resolved from window.HUB_CONFIG (set in index.html before script load).
    // Defaults to localhost so local dev works without any config.
    literatureBaseUrl: (window.HUB_CONFIG && window.HUB_CONFIG.API_BASE)
      || "http://127.0.0.1:8766",
    literatureCache: new Map(),   // key → result obj
    literaturePending: new Set(), // key currently in-flight

    chatHistory: [],              // [{role:"user"|"agent", content, actions}]
    chatPending: false,

    deZ: null,                // current disease asset (cached)
    deZCache: new Map(),      // disease → asset
    payload: null,            // current payload
    payloadCache: new Map(),  // "disease__cs" → payload

    spec: null,               // {vsAll: {padj,passes}, vsRef: {...}, geneCategory: Map}

    cy: null,
    baselineZoom: 1,
    visitedNodes: new Set(),
    pathSequence: [],
    pinnedNode: null,
    activeNode: null,
    activeStack: [],
    activeIdx: 0,
    lastCursor: { graph: null, screen: null },
  };

  // Constants matching docs/hub_payload_schema.md and concentric_egrn.py
  const R_TF = 1.8;
  const R_GENE_DELTA = 1.2;
  const POS_SCALE = 60;
  const ZOOM_SHRINK_FLOOR = 0.55;
  const ZOOM_SHRINK_POWER = 0.6;

  // ---------- DOM helpers ----------
  const $ = (id) => document.getElementById(id);
  const banner = (msg, hideAfterMs) => {
    let b = $("status-banner");
    if (!b) {
      b = document.createElement("div");
      b.id = "status-banner";
      document.body.appendChild(b);
    }
    b.textContent = msg;
    b.style.display = msg ? "block" : "none";
    if (msg && hideAfterMs) setTimeout(() => { b.style.display = "none"; }, hideAfterMs);
  };

  // ---------- fetch helpers ----------
  async function fetchJson(url) {
    const resp = await fetch(url, { cache: "force-cache" });
    if (!resp.ok) throw new Error(`fetch ${url} failed: ${resp.status}`);
    return await resp.json();
  }
  async function fetchGzJson(url) {
    const resp = await fetch(url, { cache: "force-cache" });
    if (!resp.ok) throw new Error(`fetch ${url} failed: ${resp.status}`);
    const buf = new Uint8Array(await resp.arrayBuffer());
    // fflate is loaded via UMD as `fflate` global
    const decompressed = fflate.gunzipSync(buf);
    const txt = new TextDecoder("utf-8").decode(decompressed);
    return JSON.parse(txt);
  }

  // ---------- init ----------
  async function init() {
    banner("Loading manifest…");
    try {
      S.manifest = await fetchJson("manifest.json");
      S.lineage = await fetchJson("lineage.json");
    } catch (err) {
      banner(`Manifest load failed: ${err.message}`);
      return;
    }
    populateDiseaseDropdown();
    populatePresetDropdown();
    wireUi();
    const firstDisease = Object.keys(S.manifest.diseases)[0];
    await selectDisease(firstDisease);
    banner("");
  }

  function populateDiseaseDropdown() {
    const sel = $("disease-select");
    sel.innerHTML = "";
    for (const [id, d] of Object.entries(S.manifest.diseases)) {
      const opt = document.createElement("option");
      opt.value = id; opt.textContent = d.label;
      sel.appendChild(opt);
    }
  }
  function populatePresetDropdown() {
    const sel = $("preset-select");
    while (sel.options.length > 1) sel.remove(1);
    for (const p of (S.manifest.presets || [])) {
      const opt = document.createElement("option");
      opt.value = p.label; opt.textContent = p.label;
      sel.appendChild(opt);
    }
  }

  function populateTargetCsDropdown(filter) {
    const sel = $("target-cs-select");
    const csList = S.manifest.diseases[S.diseaseId].target_cs_list;
    const q = (filter || "").trim().toLowerCase();
    const lineage = S.lineage?.groups || {};
    const nCells = S.lineage?.n_cells || {};

    // Build cs → cell_type map from lineage
    const csToType = new Map();
    for (const [grp, members] of Object.entries(lineage))
      for (const cs of members) csToType.set(cs, grp);

    sel.innerHTML = "";
    // Group cs by their parent cell_type, sorted alphabetically
    const grouped = new Map();
    for (const cs of csList) {
      if (q && !cs.toLowerCase().includes(q)) continue;
      const t = csToType.get(cs) || "Other";
      if (!grouped.has(t)) grouped.set(t, []);
      grouped.get(t).push(cs);
    }
    const types = [...grouped.keys()].sort();
    let visibleCount = 0;
    for (const t of types) {
      const og = document.createElement("optgroup");
      og.label = t;
      for (const cs of grouped.get(t).sort()) {
        const opt = document.createElement("option");
        const n = nCells[cs];
        opt.value = cs;
        opt.textContent = `${cs}${n != null ? ` (n=${n.toLocaleString()})` : ""}`;
        og.appendChild(opt);
        visibleCount++;
      }
      sel.appendChild(og);
    }
    $("target-cs-meta").textContent = q
      ? `${visibleCount} / ${csList.length} match`
      : `${csList.length} cell_states across ${types.length} cell_types`;
    if (S.targetCs && csList.includes(S.targetCs)) sel.value = S.targetCs;
  }

  function populateScoreKeyDropdown() {
    const sel = $("score-key-select");
    sel.innerHTML = "";
    const d = S.manifest.diseases[S.diseaseId];
    for (const k of d.score_keys) {
      const opt = document.createElement("option");
      opt.value = k; opt.textContent = k;
      sel.appendChild(opt);
    }
    sel.value = d.default_score_key;
    S.scoreKey = d.default_score_key;
  }

  function populateRefPicker() {
    if (!S.deZ || !S.lineage) return;
    const csList = S.deZ.cs;
    const groups = S.lineage.groups;
    const nCells = S.lineage.n_cells || {};
    const root = $("ref-checklist");
    root.innerHTML = "";
    const seen = new Set();
    for (const [grp, members] of Object.entries(groups)) {
      const inGroup = members.filter((m) => csList.includes(m));
      if (!inGroup.length) continue;
      const hdr = document.createElement("div");
      hdr.className = "ref-group";
      hdr.innerHTML = `${grp} <span class="meta">(${inGroup.length})</span>
        <a href="#" class="ref-group-toggle" data-cts="${inGroup.join('|')}">[all]</a>
        <a href="#" class="ref-group-toggle ref-group-clear" data-cts="${inGroup.join('|')}">[none]</a>`;
      root.appendChild(hdr);
      for (const cs of inGroup) {
        seen.add(cs);
        root.appendChild(refCheckRow(cs, nCells[cs]));
      }
    }
    const orphans = csList.filter((c) => !seen.has(c));
    if (orphans.length) {
      const hdr = document.createElement("div");
      hdr.className = "ref-group";
      hdr.textContent = `Other (${orphans.length})`;
      root.appendChild(hdr);
      for (const cs of orphans) root.appendChild(refCheckRow(cs, nCells[cs]));
    }
    // Bulk-select handlers
    root.querySelectorAll(".ref-group-toggle").forEach((a) => {
      a.addEventListener("click", (ev) => {
        ev.preventDefault();
        const cts = ev.target.dataset.cts.split("|");
        const clear = ev.target.classList.contains("ref-group-clear");
        for (const cs of cts) {
          if (cs === S.targetCs) continue;
          if (clear) S.refCsList = S.refCsList.filter((x) => x !== cs);
          else if (!S.refCsList.includes(cs)) S.refCsList.push(cs);
        }
        // Sync checkboxes
        root.querySelectorAll(".ref-item input").forEach((cb) => {
          cb.checked = S.refCsList.includes(cb.value);
        });
        updateRefCount();
        if (S.refMode === "selected") debouncedRecomputeSpec();
      });
    });
    updateRefCount();
  }
  function refCheckRow(cs, nCells) {
    const row = document.createElement("label");
    row.className = "ref-item";
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.value = cs;
    cb.checked = S.refCsList.includes(cs);
    if (cs === S.targetCs) {
      cb.disabled = true;
      row.classList.add("target-cs");
    }
    cb.addEventListener("change", () => {
      if (cb.checked) {
        if (!S.refCsList.includes(cs)) S.refCsList.push(cs);
      } else {
        S.refCsList = S.refCsList.filter((x) => x !== cs);
      }
      updateRefCount();
      if (S.refMode === "selected") debouncedRecomputeSpec();
    });
    const span = document.createElement("span");
    const tag = cs === S.targetCs ? " (target)" : "";
    const n = nCells != null ? ` <span class="meta">n=${nCells.toLocaleString()}</span>` : "";
    span.innerHTML = `${cs}${tag}${n}`;
    row.appendChild(cb);
    row.appendChild(span);
    return row;
  }
  function updateRefCount() {
    $("ref-count").textContent = S.refCsList.length.toString();
  }

  function populateEdgeCtSelect() {
    const sel = $("edge-ct-select");
    sel.innerHTML = "";
    if (!S.payload) return;
    const cts = new Set();
    for (const e of S.payload.edges) for (const c of e.cts) cts.add(c);
    for (const c of [...cts].sort()) {
      const opt = document.createElement("option");
      opt.value = c; opt.textContent = c;
      sel.appendChild(opt);
    }
  }

  // ---------- selection actions ----------
  async function selectDisease(id) {
    if (S.diseaseId === id) return;
    S.diseaseId = id;
    $("disease-select").value = id;
    populateScoreKeyDropdown();
    populateTargetCsDropdown();
    // Load DE-z asset (cached)
    if (!S.deZCache.has(id)) {
      const file = S.manifest.diseases[id].de_z_file;
      banner(`Fetching DE-z matrix for ${id}…`);
      console.time("fetch deZ");
      const asset = file.endsWith(".gz") ? await fetchGzJson(file) : await fetchJson(file);
      console.timeEnd("fetch deZ");

      const nGenes = asset.genes.length;
      const nCs = asset.cs.length;
      banner(`Indexing ${nGenes.toLocaleString()} × ${nCs} DE-z matrix…`);

      // Build flat Float32Array in chunks so the spinner can paint.
      console.time("build deZ flat");
      const flat = new Float32Array(nGenes * nCs);
      const CHUNK = 2000;
      let i = 0;
      while (i < nGenes) {
        const end = Math.min(i + CHUNK, nGenes);
        for (let r = i; r < end; r++) {
          const row = asset.de_z_matrix[r];
          for (let j = 0; j < nCs; j++) flat[r * nCs + j] = row[j];
        }
        i = end;
        if (i < nGenes) await new Promise((res) => setTimeout(res, 0));   // yield
      }
      console.timeEnd("build deZ flat");

      const csIdx = new Map();
      asset.cs.forEach((c, i) => csIdx.set(c, i));
      const geneIdx = new Map();
      asset.genes.forEach((g, i) => geneIdx.set(g, i));
      // Drop the nested array to free ~10 MB before stashing in cache
      delete asset.de_z_matrix;
      S.deZCache.set(id, {
        ...asset, flat, csIdx, geneIdx, nGenes, nCs,
      });
    }
    S.deZ = S.deZCache.get(id);
    const firstCs = S.manifest.diseases[id].target_cs_list[0];
    await selectTargetCs(firstCs);
  }

  async function selectTargetCs(cs) {
    S.targetCs = cs;
    $("target-cs-select").value = cs;
    const key = `${S.diseaseId}__${cs}`;
    if (!S.payloadCache.has(key)) {
      banner(`Fetching ${cs} payload…`);
      console.time(`fetch payload ${cs}`);
      const p = await fetchJson(`payloads/${key}.json`);
      console.timeEnd(`fetch payload ${cs}`);
      S.payloadCache.set(key, p);
    }
    banner("Building network…");
    S.payload = S.payloadCache.get(key);
    // Set defaults from payload only on first cs select per disease
    Object.assign(S, {
      topN: S.payload.defaults.top_n_genes,
      peaksPerSeed: S.payload.defaults.n_upstream_per_seed,
      tfsPerPeak: S.payload.defaults.n_tf_per_peak,
      fdrThresh: S.payload.defaults.fdr_thresh,
    });
    syncControlsFromState();
    populateRefPicker();
    populateEdgeCtSelect();
    rebuild();
    banner("");
  }

  function syncControlsFromState() {
    $("topn-slider").value = S.topN;
    $("topn-label").textContent = S.topN;
    $("peaks-slider").value = S.peaksPerSeed;
    $("peaks-label").textContent = S.peaksPerSeed;
    $("tfs-slider").value = S.tfsPerPeak;
    $("tfs-label").textContent = S.tfsPerPeak;
    $("fdr-input").value = S.fdrThresh;
  }

  // ---------- subgraph + layout ----------
  function buildSubgraph() {
    const p = S.payload;
    const seeds = p.seed_genes.slice(0, S.topN).map((s) => s.id);
    const seedSet = new Set(seeds);

    // Build edge cell_types lookup once per call: `${s}>${t}` → cts[]
    // Used to apply cell-type-priority sort when picking top-N peaks/TFs.
    const ctsHighlight = new Set(S.edgeHighlightCts);
    const useCtsPriority = S.cellTypePriority && ctsHighlight.size > 0;
    const edgeCts = new Map();
    if (useCtsPriority) {
      for (const e of p.edges) edgeCts.set(`${e.s}>${e.t}`, e.cts);
    }
    function sortByCtsPriority(list, edgeKey) {
      // list = [id1, id2, ...] already pre-sorted by GWAS-z (builder order).
      // Re-sort: cell-type-match first, stable wrt original index.
      if (!useCtsPriority) return list;
      const rated = list.map((id, idx) => {
        const ects = edgeCts.get(edgeKey(id)) || [];
        const m = ects.some((c) => ctsHighlight.has(c)) ? 1 : 0;
        return [id, m, idx];
      });
      rated.sort((a, b) => (b[1] - a[1]) || (a[2] - b[2]));
      return rated.map((r) => r[0]);
    }

    const peaksOrdered = [];
    const peakSet = new Set();
    for (const g of seeds) {
      const orig = p.gene_to_peaks[g] || [];
      const prioritized = sortByCtsPriority(orig, (peak) => `${peak}>${g}`);
      const lst = prioritized.slice(0, S.peaksPerSeed);
      for (const peak of lst) {
        if (!peakSet.has(peak)) { peakSet.add(peak); peaksOrdered.push(peak); }
      }
    }

    const tfSet = new Set();
    for (const peak of peakSet) {
      const orig = p.peak_to_tfs[peak] || [];
      const prioritized = sortByCtsPriority(orig, (tf) => `${tf}>${peak}`);
      const lst = prioritized.slice(0, S.tfsPerPeak);
      for (const tf of lst) tfSet.add(tf);
    }
    // remove gene-as-TF (seed_genes wins)
    for (const g of seeds) tfSet.delete(g);

    // Filter edges: only edges both endpoints inside subgraph + apply edge mode
    const edges = [];
    const cts = new Set(S.edgeHighlightCts);
    for (const e of p.edges) {
      const sIn = seedSet.has(e.s) || peakSet.has(e.s) || tfSet.has(e.s);
      const tIn = seedSet.has(e.t) || peakSet.has(e.t) || tfSet.has(e.t);
      if (!sIn || !tIn) continue;
      const hit = e.cts.some((c) => cts.has(c));
      if (S.edgeMode === "subset" && cts.size > 0 && !hit) continue;
      edges.push({ ...e, hit });
    }

    // Build node maps
    const peakMap = new Map();
    for (const peak of p.peaks) if (peakSet.has(peak.id)) peakMap.set(peak.id, peak);
    const tfMap = new Map();
    for (const tf of p.tfs) if (tfSet.has(tf.id)) tfMap.set(tf.id, tf);
    const seedMap = new Map();
    for (const s of p.seed_genes.slice(0, S.topN)) seedMap.set(s.id, s);

    return {
      seeds, seedSet, seedMap,
      peaks: peaksOrdered.filter((p) => peakMap.has(p)),
      peakMap, tfSet, tfMap, edges,
    };
  }

  function computeLayout(sg) {
    // Sort seeds by chr (numeric where possible), then by id
    function chrNum(c) {
      const u = (c || "?").toString().toUpperCase();
      if (u === "X") return 23;
      if (u === "Y") return 24;
      if (u === "M") return 25;
      const n = parseInt(u, 10);
      return Number.isFinite(n) ? n : 99;
    }
    function peakSortKey(id) {
      const m = id.match(/^chr([0-9XYMxym]+):(\d+)-/);
      if (!m) return [99, 0];
      return [chrNum(m[1]), parseInt(m[2], 10)];
    }

    const seeds = [...sg.seeds].sort((a, b) => {
      const ca = chrNum(sg.seedMap.get(a)?.chr);
      const cb = chrNum(sg.seedMap.get(b)?.chr);
      return ca - cb || a.localeCompare(b);
    });
    const peaks = [...sg.peaks].sort((a, b) => {
      const ka = peakSortKey(a), kb = peakSortKey(b);
      return ka[0] - kb[0] || ka[1] - kb[1];
    });

    // Adaptive R_PEAK based on TF label length
    let maxTfChars = 5;
    for (const t of sg.tfSet) maxTfChars = Math.max(maxTfChars, t.length);
    const charW = 0.10;  // rough data-units per character
    const maxTfRadial = maxTfChars * charW;
    const R_PEAK = Math.max(3.5, R_TF + 0.25 + maxTfRadial + 0.5);
    const R_GENE = R_PEAK + R_GENE_DELTA;

    const pos = new Map();
    const geneAngles = new Map();
    seeds.forEach((g, i) => {
      const ang = Math.PI / 2 - (2 * Math.PI * i) / Math.max(1, seeds.length);
      geneAngles.set(g, ang);
      pos.set(g, { x: R_GENE * Math.cos(ang), y: R_GENE * Math.sin(ang) });
    });
    const peakAngles = new Map();
    peaks.forEach((p, i) => {
      const ang = Math.PI / 2 - (2 * Math.PI * i) / Math.max(1, peaks.length);
      peakAngles.set(p, ang);
      pos.set(p, { x: R_PEAK * Math.cos(ang), y: R_PEAK * Math.sin(ang) });
    });

    // TF barycenter: mean angle of connected peaks
    const tfBary = new Map();
    const tfPeaks = new Map();
    for (const e of sg.edges) {
      if (e.etype !== "TF2R") continue;
      if (!sg.tfSet.has(e.s)) continue;
      if (!peakAngles.has(e.t)) continue;
      if (!tfPeaks.has(e.s)) tfPeaks.set(e.s, []);
      tfPeaks.get(e.s).push(peakAngles.get(e.t));
    }
    for (const tf of sg.tfSet) {
      const angs = tfPeaks.get(tf) || [Math.PI / 2];
      // mean of unit vectors
      let sx = 0, sy = 0;
      for (const a of angs) { sx += Math.cos(a); sy += Math.sin(a); }
      tfBary.set(tf, Math.atan2(sy / angs.length, sx / angs.length));
    }
    const tfsSorted = [...sg.tfSet].sort((a, b) => tfBary.get(a) - tfBary.get(b));
    tfsSorted.forEach((t, i) => {
      const ang = Math.PI / 2 - (2 * Math.PI * i) / Math.max(1, tfsSorted.length);
      pos.set(t, { x: R_TF * Math.cos(ang), y: R_TF * Math.sin(ang) });
    });

    return { pos, R_PEAK, R_GENE };
  }

  // ---------- specificity ----------
  function computeBothSpec(sg) {
    const dz = S.deZ;
    const tgtIdx = dz.csIdx.get(S.targetCs);
    if (tgtIdx === undefined) {
      console.warn(`target_cs ${S.targetCs} not in de_z`);
      return null;
    }
    const refAllIdx = new Int32Array(dz.nCs - 1);
    let w = 0;
    for (let j = 0; j < dz.nCs; j++) if (j !== tgtIdx) refAllIdx[w++] = j;
    const vsAll = SpecificityMath.computeSpecificity(
      dz.flat, dz.nGenes, dz.nCs, tgtIdx, refAllIdx, S.fdrThresh,
    );

    let refIdxArr;
    let refFellBack = false;
    if (S.refMode === "all") {
      refIdxArr = refAllIdx;
    } else {
      const idxs = [];
      for (const c of S.refCsList) {
        if (c === S.targetCs) continue;
        const i = dz.csIdx.get(c);
        if (i !== undefined) idxs.push(i);
      }
      if (idxs.length < 2) {
        refIdxArr = refAllIdx;
        refFellBack = true;
        banner(`Need ≥2 reference cell_states (got ${idxs.length}). Falling back to vs ALL.`, 4000);
      } else {
        refIdxArr = Int32Array.from(idxs);
      }
    }
    S.refFellBack = refFellBack;
    const vsRef = SpecificityMath.computeSpecificity(
      dz.flat, dz.nGenes, dz.nCs, tgtIdx, refIdxArr, S.fdrThresh,
    );

    // Build per-gene category map for genes in subgraph (seeds + tfs that exist in deZ)
    const cat = new Map();
    function classify(id) {
      const gi = dz.geneIdx.get(id);
      if (gi === undefined) return "neither";
      const a = vsAll.passesFdr[gi] === 1;
      const r = vsRef.passesFdr[gi] === 1;
      if (a && r) return "both";
      if (a) return "all_only";
      if (r) return "ref_only";
      return "neither";
    }
    for (const g of sg.seeds) cat.set(g, classify(g));
    for (const t of sg.tfSet) cat.set(t, classify(t));
    return {
      vsAll, vsRef, geneCategory: cat,
      nVsAll: countTrue(vsAll.passesFdr),
      nVsRef: countTrue(vsRef.passesFdr),
      dfRef: vsRef.df,
    };
  }
  function countTrue(arr) {
    let c = 0; for (let i = 0; i < arr.length; i++) if (arr[i]) c++;
    return c;
  }

  // ---------- render ----------
  function colorScaleYlOrRd(t) {
    // matplotlib YlOrRd — kept in sync with concentric_egrn.py static figures
    // so the hub viewer color encoding matches publication-ready PDFs/PNGs.
    const stops = ["#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c",
                   "#fc4e2a", "#e31a1c", "#bd0026", "#800026"];
    const x = Math.max(0, Math.min(1, t));
    const seg = x * (stops.length - 1);
    const i = Math.floor(seg);
    const f = seg - i;
    if (i >= stops.length - 1) return stops[stops.length - 1];
    return lerpHex(stops[i], stops[i + 1], f);
  }
  function lerpHex(a, b, t) {
    const ah = parseInt(a.slice(1), 16), bh = parseInt(b.slice(1), 16);
    const ar = (ah >> 16) & 255, ag = (ah >> 8) & 255, ab = ah & 255;
    const br = (bh >> 16) & 255, bg = (bh >> 8) & 255, bb = bh & 255;
    const r = Math.round(ar + (br - ar) * t);
    const g = Math.round(ag + (bg - ag) * t);
    const b2 = Math.round(ab + (bb - ab) * t);
    return `rgb(${r},${g},${b2})`;
  }
  function specBorder(cat) {
    // `both` = near-black (highest contrast, slightly larger node); single-test
    // cats keep their distinct cool hues (sky blue / emerald green).
    return { both: "#0f172a", all_only: "#0ea5e9", ref_only: "#047857",
             neither: "#a8a29e" }[cat] || "#a8a29e";
  }

  function buildElements(sg, layout) {
    const elements = [];
    const dz = S.deZ;
    // Determine GWAS-z color scale across all visible nodes (seeds + peaks + TFs).
    // Use log10(x+1) to compress the dynamic range — score columns can run into
    // thousands while most values are tiny.
    const vals = [];
    for (const g of sg.seeds) {
      const gi = dz.geneIdx.get(g);
      if (gi !== undefined) vals.push(Math.log10(Math.max(0, dz.gwas_z[gi]) + 1));
    }
    for (const id of sg.peaks) {
      const pk = sg.peakMap.get(id);
      if (pk?.gwas_z != null) vals.push(Math.log10(Math.max(0, pk.gwas_z) + 1));
    }
    for (const tfId of sg.tfSet) {
      const tf = sg.tfMap.get(tfId);
      if (tf?.gwas_z != null) vals.push(Math.log10(Math.max(0, tf.gwas_z) + 1));
    }
    const dataMax = vals.length ? Math.max(...vals) : 1;
    const dataMin = 0;
    // Allow user override via color-scale controls
    const vmax = S.colorMaxOverride != null ? S.colorMaxOverride : dataMax;
    const vmin = S.colorMinOverride != null ? S.colorMinOverride : dataMin;
    S.colorScaleAuto = { min: dataMin, max: dataMax };
    S.colorScale = { min: vmin, max: vmax };

    // Compute degree
    const degree = new Map();
    for (const e of sg.edges) {
      degree.set(e.s, (degree.get(e.s) || 0) + 1);
      degree.set(e.t, (degree.get(e.t) || 0) + 1);
    }

    function gwasColorFromVal(rawVal) {
      const v = Math.log10(Math.max(0, rawVal ?? 0) + 1);
      const t = vmax > vmin ? (v - vmin) / (vmax - vmin) : 0;
      return colorScaleYlOrRd(t);
    }
    // Radial-gradient stops: bright core → main color → slightly darker edge
    function gwasGradient(rawVal) {
      const main = gwasColorFromVal(rawVal);
      // Light center (Cividis low-end) gives the inner highlight; main is mid;
      // edge is the same main color (the gradient itself implies depth).
      return `#ffffff ${main} ${main}`;
    }
    function geneGwasColor(id) {
      const gi = dz.geneIdx.get(id);
      return gwasColorFromVal(gi !== undefined ? dz.gwas_z[gi] : 0);
    }
    function geneGwasGradient(id) {
      const gi = dz.geneIdx.get(id);
      return gwasGradient(gi !== undefined ? dz.gwas_z[gi] : 0);
    }
    function nodeSize(d, kind, cat) {
      let base;
      if (kind === "gene") base = 28;
      else if (kind === "peak") base = 12 + Math.min(d, 12) * 1.5;
      else base = 14 + Math.min(d, 14) * 1.2;
      // `both` cat = highest-priority signal → ~1.25× larger
      return cat === "both" ? base * 1.25 : base;
    }

    // Genes
    for (const g of sg.seeds) {
      const cat = S.spec?.geneCategory.get(g) || "neither";
      const p = layout.pos.get(g);
      const seed = sg.seedMap.get(g);
      elements.push({
        data: {
          id: g, label: g, kind: "gene",
          gwas_color: geneGwasColor(g),
          gradient_stops: geneGwasGradient(g),
          border: specBorder(cat),
          cat,
          sz: nodeSize(degree.get(g) || 0, "gene", cat),
          chr: seed?.chr || "?",
          gwas_z: seed?.gwas_z,
          de_z: seed?.target_de_z,
          raw: seed,
        },
        position: { x: p.x * POS_SCALE, y: -p.y * POS_SCALE },
        classes: `gene cat-${cat}`,
      });
    }
    // Peaks
    for (const id of sg.peaks) {
      const pk = sg.peakMap.get(id);
      const p = layout.pos.get(id);
      const hasSnp = pk.snp_overlap;
      const classes = ["peak", `cat-neither`];
      if (hasSnp && S.snpEmphasis) classes.push("snp-default");
      if (hasSnp && S.snpHalo) classes.push("snp");
      elements.push({
        data: {
          id, label: id, kind: "peak",
          gwas_color: gwasColorFromVal(pk.gwas_z),
          gradient_stops: gwasGradient(pk.gwas_z),
          border: hasSnp && S.snpEmphasis ? "#2563eb" : "#94a3b8",
          cat: "neither",
          sz: nodeSize(degree.get(id) || 0, "peak", "neither"),
          chr: pk.chr,
          pos: pk.pos,
          gwas_z: pk.gwas_z,
          has_snp: hasSnp ? 1 : 0,
          snp_rsids: pk.snp_rsids,
          raw: pk,
        },
        position: { x: p.x * POS_SCALE, y: -p.y * POS_SCALE },
        classes: classes.join(" "),
      });
    }
    // TFs
    for (const tfId of sg.tfSet) {
      const tf = sg.tfMap.get(tfId);
      const cat = S.spec?.geneCategory.get(tfId) || "neither";
      const p = layout.pos.get(tfId);
      elements.push({
        data: {
          id: tfId, label: tfId, kind: "tf",
          gwas_color: gwasColorFromVal(tf?.gwas_z ?? 0),
          gradient_stops: gwasGradient(tf?.gwas_z ?? 0),
          border: specBorder(cat),
          cat,
          sz: nodeSize(degree.get(tfId) || 0, "tf", cat),
          gwas_z: tf?.gwas_z,
          raw: tf,
        },
        position: { x: p.x * POS_SCALE, y: -p.y * POS_SCALE },
        classes: `tf cat-${cat}`,
      });
    }
    // Edges
    for (const e of sg.edges) {
      elements.push({
        data: {
          source: e.s, target: e.t,
          etype: e.etype, hit: e.hit ? 1 : 0,
          cts: e.cts.join(", "),
        },
        classes: `edge ${e.etype} ${e.hit ? "edge-hit" : ""}`,
      });
    }
    return { elements, degree };
  }

  // ---------- cytoscape style ----------
  function buildStyle() {
    return [
      { selector: "node", style: {
          // Solid fill — GWAS-z magnitude as Cividis color, visible across
          // the entire node face (radial gradient was washing it out).
          "background-color": "data(gwas_color)",
          "background-opacity": 0.94,
          // Universal ultra-thin slate hairline keeps SHAPE definition for
          // every node without any heavy ring. Specificity adds a stronger
          // colored ring on top via the .cat-* rules.
          "border-color": "data(border)",
          "border-width": 1,
          "border-opacity": 0.22,
          "label": "data(label)",
          "font-size": 8,
          "font-weight": 500,
          "color": "#334155",                       // text fill — slate
          "text-valign": "center",
          "text-halign": "center",
          "text-outline-color": "#f7f8fa",
          "text-outline-width": 1.4,
          "text-opacity": 0,
          "letter-spacing": 0.2,
          "width": "data(sz)",
          "height": "data(sz)",
          "transition-property": "border-width, border-color, border-opacity, underlay-opacity, font-size, opacity",
          "transition-duration": 120,
        } },
      { selector: "node.gene", style: {
          "font-size": 10, "font-weight": 600,
          "text-opacity": S.showGeneLabels ? 1 : 0, "z-index": 20,
          // Warm stone hairline — visible on both pale-yellow and dark-crimson
          // fills without ever reading as "black".
          "border-width": 1, "border-color": "#a8a29e", "border-opacity": 0.55,
        } },
      { selector: "node.tf", style: {
          "z-index": 10,
          "border-width": 1, "border-color": "#a8a29e", "border-opacity": 0.55,
        } },
      { selector: "node.peak", style: {
          "shape": "round-rectangle", "z-index": 15,
          "border-width": 0.8, "border-color": "#a8a29e", "border-opacity": 0.5,
        } },
      // Specific nodes: stronger colored ring. Genes (outer ring) get the
      // boldest treatment since they're the figure's main story-tellers.
      { selector: "node.cat-both",     style: { "border-color": "#0f172a", "border-opacity": 1, "border-width": 3.2 } },
      { selector: "node.cat-all_only", style: { "border-color": "#0ea5e9", "border-opacity": 1, "border-width": 3.2 } },
      { selector: "node.cat-ref_only", style: { "border-color": "#047857", "border-opacity": 1, "border-width": 3.2 } },
      // Gene (outer ring) gets a touch more weight; cat-both is also slightly
      // larger via the `sz` data attribute set at build time.
      { selector: "node.gene.cat-both",     style: { "border-width": 4 } },
      { selector: "node.gene.cat-all_only", style: { "border-width": 4 } },
      { selector: "node.gene.cat-ref_only", style: { "border-width": 4 } },
      { selector: "node.peak.snp-default", style: {
          "border-color": "#2563eb", "border-opacity": 0.7,
          "border-width": 1.2, "z-index": 18,
          "underlay-color": "#bfdbfe", "underlay-opacity": 0.35, "underlay-padding": 2,
        } },
      { selector: "node.snp", style: {
          "overlay-color": "#2563eb", "overlay-opacity": 0.16, "overlay-padding": 5,
        } },
      { selector: "node.no-label", style: { "text-opacity": 0 } },
      { selector: "node.force-label", style: { "text-opacity": 1 } },
      { selector: "node.no-label.force-label", style: { "text-opacity": 1, "z-index": 5 } },
      // Labelled TFs should rise above peaks (15-18) and genes (20) so the
      // text doesn't get clipped by overlapping rings.
      { selector: "node.tf.force-label", style: { "z-index": 22 } },
      { selector: "node.tf.force-label.cat-both",     style: { "z-index": 24 } },
      { selector: "node.tf.force-label.cat-all_only", style: { "z-index": 24 } },
      { selector: "node.tf.force-label.cat-ref_only", style: { "z-index": 24 } },
      { selector: "edge", style: {
          "width": 0.55,
          "line-color": "#94a3b8",
          "curve-style": "bezier",
          "control-point-step-size": 30,
          "opacity": 0.32,
          "target-arrow-shape": "none",
          "transition-property": "opacity, width, line-color",
          "transition-duration": 120,
        } },
      { selector: "edge.edge-hit", style: {
          "line-color": "#f97316", "width": 1.1, "opacity": 0.75, "z-index": 5,
        } },
      { selector: ".faded", style: { opacity: 0.07, "text-opacity": 0 } },
      { selector: "node.visited", style: {
          "underlay-color": "#fcd34d", "underlay-opacity": 0.45, "underlay-padding": 5,
          "opacity": 1, "text-opacity": 1, "z-index": 25,
        } },
      { selector: "edge.path-edge", style: {
          "line-color": "#f59e0b", "opacity": 0.92, "width": 2.4, "z-index": 30,
        } },
      { selector: "node.hovered", style: {
          "z-index": 9999,
          "border-width": 1.2, "border-color": "#f59e0b", "border-opacity": 0.92,
          "underlay-color": "#fbbf24", "underlay-opacity": 0.32, "underlay-padding": 4,
          "background-opacity": 1, "opacity": 1, "text-opacity": 1,
        } },
      { selector: "node.pinned", style: {
          "z-index": 10000,
          "border-width": 1.4, "border-color": "#2563eb", "border-opacity": 0.95,
          "underlay-color": "#60a5fa", "underlay-opacity": 0.34, "underlay-padding": 5,
          "background-opacity": 1, "opacity": 1, "text-opacity": 1,
        } },
      { selector: "node.hovered.visited", style: {
          "border-color": "#dc2626", "border-width": 1.4, "border-opacity": 0.9,
          "underlay-color": "#f87171", "underlay-opacity": 0.4, "underlay-padding": 4,
        } },
      // Chat-driven highlights — soft magenta glow, hairline border only
      { selector: "node.chat-highlight", style: {
          "border-color": "#ec4899", "border-width": 1.4, "border-opacity": 0.9,
          "underlay-color": "#f9a8d4", "underlay-opacity": 0.55, "underlay-padding": 6,
          "opacity": 1, "text-opacity": 1, "z-index": 35,
        } },
      { selector: "edge.chat-highlight", style: {
          "line-color": "#ec4899", "width": 1.8, "opacity": 0.92, "z-index": 35,
        } },
    ];
  }

  // ---------- TF label policy ----------
  function applyTfLabelMode(elements, sg) {
    // tag each TF node with 'no-label' OR 'force-label' according to mode
    const degree = new Map();
    for (const e of sg.edges) {
      if (e.etype !== "TF2R") continue;
      degree.set(e.s, (degree.get(e.s) || 0) + 1);
    }
    const tfsByDeg = [...sg.tfSet].sort((a, b) => (degree.get(b) || 0) - (degree.get(a) || 0));
    const topSet = new Set(tfsByDeg.slice(0, S.tfTopN));
    const cy = S.cy;
    cy.batch(() => {
      cy.nodes("[kind = 'tf']").forEach((n) => {
        n.removeClass("no-label force-label");
        const cat = n.data("cat");
        const id = n.id();
        const isTop = topSet.has(id);
        const isSpec = cat === "both" || cat === "all_only" || cat === "ref_only";
        let show = false;
        switch (S.tfLabelMode) {
          case "all": show = true; break;
          case "none": show = false; break;
          case "specific": show = isSpec; break;
          case "top": show = isTop; break;
          case "top_specific": default: show = isTop || isSpec; break;
        }
        if (show) n.addClass("force-label");
        else n.addClass("no-label");
      });
      cy.nodes("[kind = 'peak']").forEach((n) => {
        n.removeClass("no-label force-label");
        if (!S.showPeakLabels) n.addClass("no-label");
      });
      cy.nodes("[kind = 'gene']").forEach((n) => {
        n.removeClass("no-label force-label");
        if (!S.showGeneLabels) n.addClass("no-label");
      });
    });
  }

  // ---------- main rebuild ----------
  let _rebuildScheduled = null;
  function scheduleRebuild() {
    if (_rebuildScheduled) clearTimeout(_rebuildScheduled);
    _rebuildScheduled = setTimeout(rebuild, 80);
  }
  let _specScheduled = null;
  function debouncedRecomputeSpec() {
    if (_specScheduled) clearTimeout(_specScheduled);
    _specScheduled = setTimeout(recomputeSpecOnly, 80);
  }

  function rebuild() {
    if (!S.payload || !S.deZ) return;
    const sg = buildSubgraph();
    const layout = computeLayout(sg);
    S.spec = computeBothSpec(sg);
    const { elements } = buildElements(sg, layout);
    renderCy(elements);
    applyTfLabelMode(elements, sg);
    updateStats(sg);
    updateColorScaleUi();
  }

  function updateColorScaleUi() {
    const auto = S.colorScaleAuto, cur = S.colorScale;
    const fmt = (v) => v.toFixed(2);
    $("color-min-lbl").textContent = fmt(cur.min);
    $("color-mid-lbl").textContent = `${fmt((cur.min + cur.max) / 2)} (log10 GWAS-z+1)`;
    $("color-max-lbl").textContent = fmt(cur.max);
    if (S.colorMinOverride == null) $("color-min-input").placeholder = `auto: ${fmt(auto.min)}`;
    if (S.colorMaxOverride == null) $("color-max-input").placeholder = `auto: ${fmt(auto.max)}`;
  }

  function recomputeSpecOnly() {
    if (!S.payload || !S.deZ || !S.cy) return;
    const sg = buildSubgraph();   // cheap, same as current view
    S.spec = computeBothSpec(sg);
    // Update node cat/border in-place
    S.cy.batch(() => {
      S.cy.nodes().forEach((n) => {
        const id = n.id();
        const cat = S.spec.geneCategory.get(id);
        if (cat === undefined) return;
        const prev = n.data("cat");
        n.data("cat", cat);
        n.data("border", specBorder(cat));
        n.removeClass(`cat-${prev}`).addClass(`cat-${cat}`);
      });
    });
    applyTfLabelMode(null, sg);
    updateStats(sg);
  }

  function renderCy(elements) {
    if (S.cy) {
      S.cy.destroy();
      S.cy = null;
    }
    S.cy = cytoscape({
      container: $("cy"),
      elements,
      style: buildStyle(),
      layout: { name: "preset", animate: false, fit: true, padding: 30 },
      wheelSensitivity: 0.2,
    });
    S.baselineZoom = S.cy.zoom();
    wireCyHandlers();
    applyZoomSizing();
  }

  function applyZoomSizing() {
    if (!S.cy) return;
    const ratio = Math.max(1, S.cy.zoom() / Math.max(S.baselineZoom, 1e-6));
    const sizeScale = Math.max(ZOOM_SHRINK_FLOOR, 1 / Math.pow(ratio, ZOOM_SHRINK_POWER));
    const lin = 1 / ratio;
    S.cy.batch(() => {
      S.cy.nodes().forEach((n) => {
        const baseSize = n.data("sz");
        const baseFontSize = n.data("kind") === "gene" ? 9 : 7;
        n.style({
          width: baseSize * sizeScale,
          height: baseSize * sizeScale,
          "font-size": baseFontSize * lin,
          "text-outline-width": 0.6 * lin,
        });
      });
    });
    S.cy.style()
      .selector("node").style({ "border-width": 1 * lin })
      .selector("node.peak").style({ "border-width": 0.8 * lin })
      .selector("node.gene").style({ "border-width": 1.5 * lin })
      .selector("node.tf").style({ "border-width": 1 * lin })
      .selector("node.cat-both").style({ "border-width": 3.2 * lin })
      .selector("node.cat-all_only").style({ "border-width": 3.2 * lin })
      .selector("node.cat-ref_only").style({ "border-width": 3.2 * lin })
      .selector("node.gene.cat-both").style({ "border-width": 4 * lin })
      .selector("node.gene.cat-all_only").style({ "border-width": 4 * lin })
      .selector("node.gene.cat-ref_only").style({ "border-width": 4 * lin })
      .selector("node.peak.snp-default").style({ "border-width": 1.2 * lin })
      .selector("node.hovered").style({ "border-width": 1.8 * lin })
      .selector("node.pinned").style({ "border-width": 2.2 * lin })
      .selector("node.hovered.visited").style({ "border-width": 2 * lin })
      .selector("node.chat-highlight").style({ "border-width": 2 * lin })
      .selector("edge").style({ "width": 0.7 * lin })
      .selector("edge.edge-hit").style({ "width": 1.3 * lin })
      .selector("edge.path-edge").style({ "width": 2.4 * lin })
      .selector("edge.chat-highlight").style({ "width": 2.2 * lin })
      .update();
  }

  // ---------- interactive handlers (hover/pin/Tab/path) ----------
  function wireCyHandlers() {
    const cy = S.cy;
    cy.on("zoom", applyZoomSizing);

    cy.on("mousemove", (evt) => {
      const pos = evt.position || evt.cyPosition;
      S.lastCursor.graph = pos;
      const cands = collectCandidates(pos);
      S.activeStack = cands;
      S.activeIdx = 0;
      setActiveNode(cands[0] || null);
    });

    cy.on("tap", "node", (evt) => {
      const n = evt.target;
      toggleVisited(n);
    });
    cy.on("tap", (evt) => {
      if (evt.target === cy) {
        S.visitedNodes.clear();
        S.pathSequence = [];
        setPinnedNode(null);
        applyVisitedContext();
      }
    });
    // keydown is attached once in wireUi (not here — wireCyHandlers runs on
    // every rebuild and would leak listeners).
  }

  function collectCandidates(cursorPos) {
    if (!cursorPos || !S.cy) return [];
    const cy = S.cy;
    const radius = 24 / cy.zoom();   // graph-coord radius scaled by zoom
    const nodes = cy.nodes().toArray();
    const inRange = [];
    for (const n of nodes) {
      const p = n.position();
      const dx = p.x - cursorPos.x;
      const dy = p.y - cursorPos.y;
      const d2 = dx * dx + dy * dy;
      if (d2 <= radius * radius) inRange.push([d2, n]);
    }
    inRange.sort((a, b) => a[0] - b[0]);
    return inRange.map((x) => x[1]);
  }

  function setActiveNode(n) {
    if (!S.cy) return;
    S.cy.nodes(".hovered").removeClass("hovered");
    if (n) n.addClass("hovered");
    S.activeNode = n;
  }

  function setPinnedNode(n) {
    if (!S.cy) return;
    S.cy.nodes(".pinned").removeClass("pinned");
    if (n) {
      n.addClass("pinned");
      showDetailPanel(n);
    } else {
      hideDetailPanel();
    }
    S.pinnedNode = n;
  }

  function toggleVisited(n) {
    const id = n.id();
    if (S.visitedNodes.has(id)) {
      // un-visit
      S.visitedNodes.delete(id);
      S.pathSequence = S.pathSequence.filter((x) => x !== id);
      if (S.pinnedNode && S.pinnedNode.id() === id) setPinnedNode(null);
    } else {
      S.visitedNodes.add(id);
      S.pathSequence.push(id);
      setPinnedNode(n);
    }
    applyVisitedContext();
  }

  function applyVisitedContext() {
    if (!S.cy) return;
    const cy = S.cy;
    if (!S.visitedNodes.size) {
      cy.batch(() => {
        cy.elements().removeClass("faded visited path-edge");
      });
      return;
    }
    const context = new Set(S.visitedNodes);
    for (const id of S.visitedNodes) {
      const n = cy.getElementById(id);
      if (!n.length) continue;
      n.neighborhood("node").forEach((nb) => context.add(nb.id()));
    }
    cy.batch(() => {
      cy.nodes().removeClass("visited faded");
      cy.edges().removeClass("path-edge faded");
      cy.nodes().forEach((n) => {
        if (S.visitedNodes.has(n.id())) n.addClass("visited");
        else if (!context.has(n.id())) n.addClass("faded");
      });
      cy.edges().forEach((e) => {
        const inCtx = context.has(e.source().id()) && context.has(e.target().id());
        if (!inCtx) e.addClass("faded");
      });
      // Path edges connect consecutive clicks
      for (let i = 1; i < S.pathSequence.length; i++) {
        const a = S.pathSequence[i - 1], b = S.pathSequence[i];
        cy.edges().forEach((e) => {
          const ids = [e.source().id(), e.target().id()];
          if (ids.includes(a) && ids.includes(b)) e.addClass("path-edge").removeClass("faded");
        });
      }
    });
  }

  // ---------- detail panel ----------
  let _detailPanel = null;
  function ensureDetailPanel() {
    if (_detailPanel) return _detailPanel;
    const panel = document.createElement("div");
    panel.className = "detail-panel";
    panel.innerHTML = `<div class="detail-header"><span id="dp-title"></span><span class="close-btn">×</span></div>
                       <div class="detail-body" id="dp-body"></div>`;
    document.body.appendChild(panel);
    panel.querySelector(".close-btn").addEventListener("click", () => {
      panel.style.display = "none";
      // Closing the panel just hides the annotation; trail stays.
    });
    const header = panel.querySelector(".detail-header");
    let drag = null;
    header.addEventListener("mousedown", (e) => {
      const r = panel.getBoundingClientRect();
      drag = { dx: e.clientX - r.left, dy: e.clientY - r.top };
      panel.classList.add("dragging");
    });
    document.addEventListener("mousemove", (e) => {
      if (!drag) return;
      panel.style.left = `${e.clientX - drag.dx}px`;
      panel.style.top = `${e.clientY - drag.dy}px`;
      panel.style.right = "auto";
    });
    document.addEventListener("mouseup", () => { drag = null; panel.classList.remove("dragging"); });
    _detailPanel = panel;
    return panel;
  }
  function showDetailPanel(n) {
    const panel = ensureDetailPanel();
    panel.style.display = "block";
    panel.querySelector("#dp-title").textContent = n.id();
    panel.querySelector("#dp-body").innerHTML = buildDetailHtml(n);
    // Wire "Ask agent" shortcuts (per-node and per-SNP)
    panel.querySelectorAll(".ask-agent-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        const kind = btn.dataset.kind;
        const key = btn.dataset.key;
        openChatAndPrefill(kind, key);
      });
    });
    // Wire AG1 variant-effect buttons
    panel.querySelectorAll(".ag1-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        openAg1Panel(btn.dataset.variant, btn.dataset.rsid);
      });
    });
    // Wire connection ID clicks → pin that node + scroll/center on it
    panel.querySelectorAll(".conn-id").forEach((span) => {
      span.addEventListener("click", () => {
        const id = span.dataset.id;
        const target = S.cy?.getElementById(id);
        if (target && target.length) {
          toggleVisited(target);
          S.cy.animate({ center: { eles: target }, zoom: Math.max(S.cy.zoom(), 0.6) },
                       { duration: 250 });
        }
      });
    });
  }
  function hideDetailPanel() {
    if (_detailPanel) _detailPanel.style.display = "none";
  }

  function buildDetailHtml(n) {
    const d = n.data();
    const kind = d.kind;
    const cat = d.cat || "neither";
    let html = `<div class="tt-sub">${kind.toUpperCase()} • <span class="tt-tag tt-tag-${cat}">${cat}</span></div>`;
    if (d.gwas_z !== undefined) html += `<div class="tt-row"><b>GWAS-z</b> ${fmt(d.gwas_z)}</div>`;
    if (d.de_z !== undefined) html += `<div class="tt-row"><b>${S.targetCs} DE-z</b> ${fmt(d.de_z)}</div>`;
    if (d.chr) html += `<div class="tt-row"><b>chr</b> ${d.chr}${d.pos ? ":" + d.pos : ""}</div>`;
    if (kind === "peak" && d.has_snp) {
      html += `<div class="tt-snps"><b>SNP-overlapping peak</b>`;
      const snps = (S.payload.snps || []).filter((s) => s.peak === d.id);
      for (const s of snps) {
        html += renderSnpBlock(s);
      }
      html += `</div>`;
    }
    if (kind === "gene" || kind === "tf") {
      html += renderGeneLinks(d.id);
    }
    // Network connections (direct + 2-hop)
    html += renderNetworkConnections(n);
    // "Ask agent →" shortcut — pre-fills the chat with a context query.
    html += renderAskAgentShortcut(kind, d);
    return html;
  }

  function renderNetworkConnections(n) {
    if (!n || !n.length) return "";
    const d = n.data();
    const kind = d.kind;
    const inc = n.incomers("node");
    const out = n.outgoers("node");

    // Build summary of a node list (for compact display)
    const fmtIdList = (arr, max = 12) => {
      const ids = arr.map((x) => x.id());
      const head = ids.slice(0, max).map((id) =>
        `<span class="conn-id" data-id="${id}">${id}</span>`).join(", ");
      const more = ids.length > max ? ` <span class="meta">+${ids.length - max} more</span>` : "";
      return head + more;
    };

    let h = `<div class="tt-conn"><div class="tt-row"><b>Network connections</b></div>`;

    if (kind === "gene") {
      // R2G upstream peaks (direct incomers) + TF2R 2-hop TFs
      const peaks = inc.filter("[kind = 'peak']");
      const tf2hop = peaks.incomers("node[kind = 'tf']").union(peaks.incomers("node[kind = 'gene']"));
      const snpPeaks = peaks.filter("[has_snp = 1]");
      h += `<div class="tt-row"><b>Upstream peaks</b> ${peaks.length}`
        + (snpPeaks.length ? ` <span class="meta">(${snpPeaks.length} with SNP)</span>` : "") + `</div>`;
      if (peaks.length) h += `<div class="tt-row">${fmtIdList(peaks)}</div>`;
      h += `<div class="tt-row"><b>Regulating TFs (2-hop)</b> ${tf2hop.length}</div>`;
      if (tf2hop.length) h += `<div class="tt-row">${fmtIdList(tf2hop)}</div>`;
    } else if (kind === "tf") {
      // TF2R downstream peaks (direct outgoers) + R2G 2-hop genes
      const peaks = out.filter("[kind = 'peak']");
      const genes = peaks.outgoers("node[kind = 'gene']");
      h += `<div class="tt-row"><b>Bound peaks</b> ${peaks.length}</div>`;
      if (peaks.length) h += `<div class="tt-row">${fmtIdList(peaks)}</div>`;
      h += `<div class="tt-row"><b>Target genes (2-hop)</b> ${genes.length}</div>`;
      if (genes.length) h += `<div class="tt-row">${fmtIdList(genes)}</div>`;
    } else if (kind === "peak") {
      // Incoming TFs + outgoing genes
      const tfs = inc.filter("[kind = 'tf'], [kind = 'gene']");
      const genes = out.filter("[kind = 'gene']");
      h += `<div class="tt-row"><b>Bound by TFs</b> ${tfs.length}</div>`;
      if (tfs.length) h += `<div class="tt-row">${fmtIdList(tfs)}</div>`;
      h += `<div class="tt-row"><b>Target genes</b> ${genes.length}</div>`;
      if (genes.length) h += `<div class="tt-row">${fmtIdList(genes)}</div>`;
    }
    h += `</div>`;
    return h;
  }

  function renderAskAgentShortcut(kind, d) {
    const id = (kind === "peak" && d.has_snp) ? null : d.id;   // peak handled below in SNP block
    if (!id) return "";
    return `<div class="tt-lit">
      <button class="ask-agent-btn" data-kind="${kind}" data-key="${id}">💬 Ask agent about ${id}</button>
    </div>`;
  }

  function renderLitResult(res) {
    const pm = (res.pubmed && res.pubmed.citations) || [];
    const epmc = (res.europe_pmc && res.europe_pmc.citations) || [];
    const ag = res.agent_citations || [];
    let h = "";
    if (res.mode) {
      h += `<div class="meta">Mode: <code>${res.mode}</code>`
        + (res.query_label ? ` · Query <code>${res.query_label}</code> (${pm.length}/${res.pubmed?.count_total || pm.length} PubMed)` : "")
        + `</div>`;
    }
    if (res.synthesis) {
      h += `<div class="lit-synthesis">${res.synthesis.replace(/\n/g, "<br>")}</div>`;
    }
    if (res.agent_summary) {
      h += `<div class="lit-synthesis"><b>Agent:</b> ${res.agent_summary.replace(/\n/g, "<br>")}</div>`;
    }
    if (res.agent_broader) {
      h += `<div class="meta"><b>Broader context:</b> ${res.agent_broader}</div>`;
    }
    if (res.agent_error) {
      h += `<div class="meta" style="color:#c62828">Agent error: ${res.agent_error}</div>`;
    }

    function renderCite(c, sourceLabel) {
      const url = c.url || (c.pmid ? `https://pubmed.ncbi.nlm.nih.gov/${c.pmid}/` : (c.doi ? `https://doi.org/${c.doi}` : ""));
      const idTxt = c.pmid ? `PMID ${c.pmid}` : (c.doi ? `DOI ${c.doi}` : "");
      const link = url ? `<a target="_blank" href="${url}">${idTxt}</a>` : idTxt;
      const sourceBadge = sourceLabel ? `<span class="lit-source">${sourceLabel}</span>` : "";
      const verifyBadge = (c.verified === false)
        ? `<span class="lit-source" style="background:#c62828">unverified</span>` : "";
      const findingLine = c.key_finding ? `<div>${c.key_finding}</div>` : "";
      return `<div class="lit-cite">${sourceBadge}${verifyBadge}`
        + `<b>${c.title || "(no title)"}</b> `
        + `<span class="meta">${c.authors || ""} <i>${c.journal || ""}</i> ${c.year || ""} ${link}</span>`
        + findingLine + `</div>`;
    }

    if (ag.length) {
      const nV = ag.filter((c) => c.verified).length;
      h += `<details class="lit-citations" open><summary>Agent citations (${nV}/${ag.length} verified)</summary>`;
      for (const c of ag) h += renderCite(c, "agent");
      h += `</details>`;
    }
    if (pm.length) {
      h += `<details class="lit-citations" ${ag.length ? "" : "open"}><summary>PubMed (${pm.length})</summary>`;
      for (const c of pm) h += renderCite(c, null);
      h += `</details>`;
    } else if (!ag.length) {
      h += `<div class="meta">No PubMed hits for this query.</div>`;
    }
    if (epmc.length) {
      h += `<details class="lit-citations"><summary>Europe PMC (${epmc.length} extra)</summary>`;
      for (const c of epmc) h += renderCite(c, c.source);
      h += `</details>`;
    }
    if (res.agent_raw) {
      h += `<details><summary>Agent raw (parse failed)</summary>`
        + `<pre style="font-size:10px;white-space:pre-wrap">${res.agent_raw.slice(0, 1000)}</pre></details>`;
    }
    return h;
  }

  async function searchLiterature(type, key, mode = "pubmed") {
    // Cache key intentionally omits mode — re-running with a different mode
    // overwrites the displayed result (one current view per node).
    const cacheKey = `${type}__${S.diseaseId}__${S.targetCs}__${key}`;
    if (S.literaturePending.has(cacheKey)) return;
    S.literaturePending.add(cacheKey);
    S.literatureCurrentMode = mode;
    if (S.pinnedNode) showDetailPanel(S.pinnedNode);

    const ctx = { disease: S.diseaseId, target_cs: S.targetCs };
    if (type === "snp") {
      const s = (S.payload.snps || []).find((x) =>
        x.rsid === key
        || (x.annotations?.ot_variant?.rsIds || []).includes(key));
      if (s) { ctx.rsid = key; ctx.chr = s.chr; ctx.pos = s.pos; }
    } else {
      ctx.symbol = key;
    }
    try {
      const resp = await fetch(`${S.literatureBaseUrl}/literature/${type}/${encodeURIComponent(key)}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ context: ctx, mode, force: true }),
      });
      if (!resp.ok) {
        const errTxt = await resp.text();
        banner(`Literature error: ${resp.status} ${errTxt.slice(0, 80)}`, 6000);
        return;
      }
      const data = await resp.json();
      S.literatureCache.set(cacheKey, data);
    } catch (err) {
      banner(`Literature backend unreachable. Start it with: python scripts/snp2cell_hub_backend.py`, 8000);
    } finally {
      S.literaturePending.delete(cacheKey);
      if (S.pinnedNode) showDetailPanel(S.pinnedNode);
    }
  }

  async function maybePreloadLiterature(type, key) {
    // GET from cache before triggering a fresh search
    const cacheKey = `${type}__${S.diseaseId}__${S.targetCs}__${key}`;
    if (S.literatureCache.has(cacheKey)) return;
    try {
      const resp = await fetch(`${S.literatureBaseUrl}/literature/${type}/${encodeURIComponent(key)}`);
      if (resp.ok) {
        S.literatureCache.set(cacheKey, await resp.json());
        if (S.pinnedNode) showDetailPanel(S.pinnedNode);
      }
    } catch (_) { /* backend not running — fine */ }
  }

  function fmt(v) { return v == null ? "—" : (Math.abs(v) > 100 ? v.toExponential(2) : Number(v).toFixed(2)); }

  function renderSnpBlock(s) {
    // s is a payload SNP entry; may include s.annotations (OT/VEP/GTEx pre-baked)
    const ann = s.annotations || {};
    const ot = ann.ot_variant || {};
    // Prefer real rs# from OT (s.rsid may be a chrN:POS fallback when GWAS
    // hitcheck has no curated rsid).
    const otRs = (ot.rsIds && ot.rsIds.length) ? ot.rsIds[0] : null;
    const displayRsid = (s.rsid && /^rs\d+/.test(s.rsid)) ? s.rsid : otRs;
    const headerId = displayRsid || s.rsid;   // fall back to chr:pos label
    let h = `<div class="tt-snp-card">`;
    h += `<div class="tt-row"><b>${headerId}</b> ${s.ref}/${s.alt} `
       + `<span class="meta">chr${s.chr}:${s.pos.toLocaleString()}</span></div>`;
    h += `<div class="tt-row"><span class="meta">PIP=${fmt(s.pip)} · `
       + `gnomAD AF=${fmt(s.gnomad_af)} · trait=${s.trait || "—"}</span></div>`;
    if (ot.mostSevereConsequence) {
      h += `<div class="tt-row"><b>Consequence</b> ${ot.mostSevereConsequence}</div>`;
    } else if (ann.most_severe_consequence) {
      h += `<div class="tt-row"><b>VEP</b> ${ann.most_severe_consequence}</div>`;
    }
    // Top transcript consequences with impact
    const tcs = (ot.transcriptConsequences || []).filter((t) => t.gene_symbol).slice(0, 4);
    if (tcs.length) {
      h += `<div class="tt-row"><b>Genes (VEP)</b> ` +
        tcs.map((t) => `<span title="${t.transcriptId || ''}">${t.gene_symbol}${t.canonical ? '*' : ''} <span class="meta">(${t.impact})</span></span>`).join(", ") +
        `</div>`;
    }
    // OT credible sets + L2G
    const cs = (ot.credibleSets || []);
    if (cs.length) {
      h += `<div class="tt-row"><b>OT credible sets</b> (${cs.length})</div>`;
      for (const c of cs.slice(0, 5)) {
        const l2g = (c.l2g_top || []).map((g) => `${g.gene_symbol} <span class="meta">L2G=${fmt(g.l2g_score)}</span>`).join(", ");
        h += `<div class="tt-row tt-cs">` +
          `${c.trait || c.study_id} <span class="meta">${c.first_author || ""}</span>` +
          (l2g ? `<br>→ ${l2g}` : "") +
          `</div>`;
      }
    }
    // GTEx heart eQTLs
    const gtex = ann.gtex_eqtl || [];
    if (gtex.length) {
      h += `<div class="tt-row"><b>GTEx heart eQTL</b><br>`
        + gtex.slice(0, 6).map((q) => `${q.gene_symbol} <span class="meta">${q.tissue.replace('Heart_', '')} p=${fmt(q.pval)} NES=${fmt(q.nes)}</span>`).join("<br>")
        + `</div>`;
    }
    h += renderSnpLinks(s, displayRsid);
    // Per-SNP "Ask agent" shortcut
    const askKey = displayRsid || s.rsid;
    h += `<div class="tt-lit">`
      + `<button class="ask-agent-btn" data-kind="snp" data-key="${askKey}">💬 Ask agent about ${askKey}</button>`
      + `</div>`;
    // AG1 variant-effect predictor (separate row — slow, opt-in)
    const variantId = `${s.chr}_${s.pos}_${s.ref}_${s.alt}`;
    h += `<div class="tt-lit">`
      + `<button class="ag1-btn" data-variant="${variantId}" data-rsid="${askKey}">🧬 AG1 variant effect (${variantId})</button>`
      + `</div>`;
    h += `</div>`;
    return h;
  }

  function renderSnpLinks(s, rsidResolved) {
    const t = S.manifest.snp_url_templates;
    const rs = rsidResolved;            // may be null if no real rs#
    function fill(tpl, useRs) {
      return tpl.replace("{rsid}", useRs || "")
                .replace("{chr}", s.chr)
                .replace("{pos}", s.pos)
                .replace("{ref}", s.ref)
                .replace("{alt}", s.alt);
    }
    const links = [];
    if (rs) {
      links.push(`<a target="_blank" href="${fill(t.dbsnp, rs)}">dbSNP</a>`);
      links.push(`<a target="_blank" href="${fill(t.gwas_catalog, rs)}">GWAS Catalog</a>`);
      links.push(`<a target="_blank" href="${fill(t.ensembl, rs)}">Ensembl</a>`);
    }
    // Coordinate-based links work for any variant
    links.push(`<a target="_blank" href="${fill(t.open_targets, rs)}">Open Targets</a>`);
    links.push(`<a target="_blank" href="${fill(t.ucsc, rs)}">UCSC</a>`);
    links.push(`<a target="_blank" href="${fill(t.gnomad, rs)}">gnomAD</a>`);
    return `<div class="tt-links">${links.join("")}</div>`;
  }

  function renderGeneLinks(symbol) {
    return `<div class="tt-links">` +
      `<a target="_blank" href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${symbol}">GeneCards</a>` +
      `<a target="_blank" href="https://www.ncbi.nlm.nih.gov/gene/?term=${symbol}">NCBI Gene</a>` +
      `<a target="_blank" href="https://www.ensembl.org/Multi/Search/Results?species=Human;q=${symbol}">Ensembl</a>` +
      `<a target="_blank" href="https://genetics.opentargets.org/?q=${symbol}">Open Targets</a>` +
      `</div>`;
  }

  // ---------- stats ----------
  function updateStats(sg) {
    const out = $("stats");
    let refLabel;
    if (S.refMode === "all") {
      refLabel = `all other (n=${S.deZ ? S.deZ.nCs - 1 : 0})`;
    } else if (S.refFellBack) {
      refLabel = `<span style="color:#c62828">selected too small → fallback to ALL</span>`;
    } else {
      refLabel = `selected (n=${S.refCsList.filter((c) => c !== S.targetCs).length})`;
    }
    out.innerHTML = `
      <div><b>Target</b> ${S.targetCs}</div>
      <div><b>Disease</b> ${S.manifest.diseases[S.diseaseId].label}</div>
      <div><b>Subgraph</b> seeds=${sg.seeds.length} peaks=${sg.peaks.length} TFs=${sg.tfSet.size} edges=${sg.edges.length}</div>
      <div><b>Spec contrast</b> vs ${refLabel} (df=${S.spec?.dfRef ?? "—"})</div>
      <div><b>FDR&lt;${S.fdrThresh}</b> vs ALL: ${S.spec?.nVsAll ?? "—"} vs ref: ${S.spec?.nVsRef ?? "—"}</div>
      <div><b>SNP peaks</b> ${sg.peaks.filter((p) => sg.peakMap.get(p)?.snp_overlap).length}</div>
    `;
  }

  // ---------- UI wiring ----------
  function wireUi() {
    // Global key handler attached once (avoid leaking on rebuilds)
    document.addEventListener("keydown", (e) => {
      if (e.key === "Tab" && S.activeStack.length > 1) {
        e.preventDefault();
        S.activeIdx = (S.activeIdx + 1) % S.activeStack.length;
        setActiveNode(S.activeStack[S.activeIdx]);
      }
    });

    $("disease-select").addEventListener("change", (e) => selectDisease(e.target.value));
    $("target-cs-select").addEventListener("change", (e) => selectTargetCs(e.target.value));
    $("target-cs-search").addEventListener("input", (e) =>
      populateTargetCsDropdown(e.target.value));
    $("score-key-select").addEventListener("change", (e) => {
      S.scoreKey = e.target.value;
      // Score key change doesn't affect DE-z, only flavour of GWAS coloring.
      // For now we ignore (single key per payload). Future: re-fetch payload.
    });

    $("preset-select").addEventListener("change", async (e) => {
      const label = e.target.value;
      if (!label) return;
      const p = (S.manifest.presets || []).find((x) => x.label === label);
      if (!p) return;
      if (p.disease !== S.diseaseId) await selectDisease(p.disease);
      if (p.target_cs !== S.targetCs) await selectTargetCs(p.target_cs);
      // Apply preset filter values
      S.topN = p.top_n_genes;
      S.edgeMode = p.edge_highlight_mode;
      S.edgeHighlightCts = p.edge_highlight_cell_types || [];
      if (p.ref_cs_list && p.ref_cs_list.length) {
        S.refMode = "selected";
        S.refCsList = [...p.ref_cs_list];
      } else {
        S.refMode = "all";
      }
      $("ref-mode").value = S.refMode;
      $("ref-picker").style.display = S.refMode === "selected" ? "block" : "none";
      $("edge-mode").value = S.edgeMode;
      populateRefPicker();
      // Sync edge ct select
      const csel = $("edge-ct-select");
      for (const opt of csel.options) opt.selected = S.edgeHighlightCts.includes(opt.value);
      syncControlsFromState();
      rebuild();
    });

    $("topn-slider").addEventListener("input", (e) => {
      S.topN = parseInt(e.target.value);
      $("topn-label").textContent = S.topN;
      scheduleRebuild();
    });
    $("peaks-slider").addEventListener("input", (e) => {
      S.peaksPerSeed = parseInt(e.target.value);
      $("peaks-label").textContent = S.peaksPerSeed;
      scheduleRebuild();
    });
    $("tfs-slider").addEventListener("input", (e) => {
      S.tfsPerPeak = parseInt(e.target.value);
      $("tfs-label").textContent = S.tfsPerPeak;
      scheduleRebuild();
    });

    $("ref-mode").addEventListener("change", (e) => {
      S.refMode = e.target.value;
      $("ref-picker").style.display = S.refMode === "selected" ? "block" : "none";
      debouncedRecomputeSpec();
    });

    $("ref-search").addEventListener("input", (e) => {
      const q = e.target.value.toLowerCase();
      document.querySelectorAll("#ref-checklist .ref-item").forEach((row) => {
        const txt = row.querySelector("span").textContent.toLowerCase();
        row.style.display = txt.includes(q) ? "" : "none";
      });
    });

    $("fdr-input").addEventListener("change", (e) => {
      S.fdrThresh = parseFloat(e.target.value);
      debouncedRecomputeSpec();
    });

    $("edge-mode").addEventListener("change", (e) => {
      S.edgeMode = e.target.value;
      scheduleRebuild();
    });
    $("edge-ct-select").addEventListener("change", (e) => {
      S.edgeHighlightCts = Array.from(e.target.selectedOptions).map((o) => o.value);
      scheduleRebuild();
    });
    $("celltype-priority-toggle").addEventListener("change", (e) => {
      S.cellTypePriority = e.target.checked;
      scheduleRebuild();
    });

    $("tf-label-mode").addEventListener("change", (e) => {
      S.tfLabelMode = e.target.value;
      if (S.cy) applyTfLabelMode(null, buildSubgraph());
    });
    $("tf-top-n").addEventListener("change", (e) => {
      S.tfTopN = parseInt(e.target.value) || 0;
      if (S.cy) applyTfLabelMode(null, buildSubgraph());
    });
    $("gene-label-toggle").addEventListener("change", (e) => {
      S.showGeneLabels = e.target.checked;
      if (S.cy) applyTfLabelMode(null, buildSubgraph());
    });
    $("peak-label-toggle").addEventListener("change", (e) => {
      S.showPeakLabels = e.target.checked;
      if (S.cy) applyTfLabelMode(null, buildSubgraph());
    });
    $("snp-emph-toggle").addEventListener("change", (e) => {
      S.snpEmphasis = e.target.checked;
      scheduleRebuild();
    });
    $("snp-halo-toggle").addEventListener("change", (e) => {
      S.snpHalo = e.target.checked;
      scheduleRebuild();
    });

    $("reset-btn").addEventListener("click", () => {
      if (S.cy) S.cy.fit(null, 40);
      S.visitedNodes.clear();
      S.pathSequence = [];
      setPinnedNode(null);
      applyVisitedContext();
    });
    $("rebuild-btn").addEventListener("click", () => rebuild());
    $("color-min-input").addEventListener("change", (e) => {
      const v = e.target.value;
      S.colorMinOverride = v === "" ? null : parseFloat(v);
      scheduleRebuild();
    });
    $("color-max-input").addEventListener("change", (e) => {
      const v = e.target.value;
      S.colorMaxOverride = v === "" ? null : parseFloat(v);
      scheduleRebuild();
    });
    $("color-reset-btn").addEventListener("click", () => {
      S.colorMinOverride = null;
      S.colorMaxOverride = null;
      $("color-min-input").value = "";
      $("color-max-input").value = "";
      scheduleRebuild();
    });

    $("export-png-btn").addEventListener("click", () => {
      if (!S.cy) return;
      const png = S.cy.png({ scale: 3, full: true, bg: "#fafafa" });
      const a = document.createElement("a");
      a.href = png;
      a.download = `snp2cell_${S.diseaseId}_${S.targetCs}.png`;
      a.click();
    });
  }

  // ---------- chat panel / Ask Agent ----------
  function buildGraphState() {
    if (!S.payload || !S.deZ) return {};
    const sg = buildSubgraph();
    const dz = S.deZ;

    // ALL visible seeds (no truncation)
    const seedRecords = sg.seeds.map((g) => {
      const seed = sg.seedMap.get(g);
      return {
        id: g, chr: seed?.chr,
        gwas_z: seed?.gwas_z,
        target_de_z: seed?.target_de_z,
        cat: S.spec?.geneCategory.get(g) || "neither",
      };
    });

    // ALL visible peaks
    const peakRecords = sg.peaks.map((id) => {
      const p = sg.peakMap.get(id);
      return {
        id, chr: p?.chr, pos: p?.pos, gwas_z: p?.gwas_z,
        snp_overlap: p?.snp_overlap ? 1 : 0,
        snp_rsids: (p?.snp_rsids || []).slice(0, 4),
      };
    });

    // ALL visible TFs
    const tfRecords = [...sg.tfSet].map((tfId) => {
      const tf = sg.tfMap.get(tfId);
      return {
        id: tfId, gwas_z: tf?.gwas_z,
        cat: S.spec?.geneCategory.get(tfId) || "neither",
      };
    });

    // ALL SNPs in payload that overlap visible peaks
    const visiblePeakIds = new Set(sg.peaks);
    const snpRecords = [];
    for (const s of S.payload.snps || []) {
      if (!visiblePeakIds.has(s.peak)) continue;
      const ann = s.annotations || {};
      const ot = ann.ot_variant || {};
      const cs = ot.credibleSets || [];
      const topL2g = cs[0]?.l2g_top || [];
      const gtex = ann.gtex_eqtl || [];
      snpRecords.push({
        rsid: (ot.rsIds && ot.rsIds[0]) || s.rsid,
        peak: s.peak,
        chr: s.chr, pos: s.pos, ref: s.ref, alt: s.alt,
        pip: s.pip, trait: s.trait,
        consequence: ot.mostSevereConsequence,
        l2g_top_genes: topL2g.slice(0, 3).map((g) => `${g.gene_symbol}(${(g.l2g_score||0).toFixed(2)})`),
        gtex_heart_eqtl_genes: [...new Set(gtex.map((q) => q.gene_symbol))].slice(0, 5),
      });
    }

    // ALL visible edges, compressed: cell_types lookup table to keep size down
    const ctTable = [];
    const ctToIdx = new Map();
    const getCtIdx = (c) => {
      if (!ctToIdx.has(c)) { ctToIdx.set(c, ctTable.length); ctTable.push(c); }
      return ctToIdx.get(c);
    };
    const edgeRecords = sg.edges.map((e) => ({
      s: e.s, t: e.t, k: e.etype === "TF2R" ? 1 : 0,
      cts: e.cts.map(getCtIdx),
    }));

    // cell_type → state hierarchy from lineage.json for context
    const lineage = S.lineage?.groups || {};
    const targetCellType = Object.entries(lineage).find(
      ([_g, members]) => members.includes(S.targetCs))?.[0] || null;

    const availableTargetCsList = S.manifest?.diseases[S.diseaseId]?.target_cs_list || [];
    const allCsInDeZ = S.deZ?.cs || [];
    const allDiseases = S.manifest ? Object.entries(S.manifest.diseases).map(
      ([id, d]) => ({ id, label: d.label, efo_id: d.efo_id, n_cs: d.target_cs_list.length })) : [];

    return {
      disease: S.diseaseId,
      target_cs: S.targetCs,
      target_cell_type: targetCellType,    // parent cell_type from h5ad
      manifest_disease_list: allDiseases,                  // valid set_disease values
      manifest_target_cs_list: availableTargetCsList,   // valid set_target_cs values
      all_cs_in_de_z: allCsInDeZ,                       // valid set_ref_cs values
      ref_mode: S.refMode,
      ref_cs_list: S.refCsList,
      edge_highlight_cell_types: S.edgeHighlightCts,
      edge_highlight_mode: S.edgeMode,
      edge_celltype_priority_on: S.cellTypePriority,
      top_n_current: S.topN,
      peaks_per_seed_current: S.peaksPerSeed,
      tfs_per_peak_current: S.tfsPerPeak,
      fdr_thresh: S.fdrThresh,
      spec_passes_vs_all: S.spec?.nVsAll,
      spec_passes_vs_ref: S.spec?.nVsRef,
      cell_type_table: ctTable,             // index → name (referenced by edges[].cts)
      n_seeds: seedRecords.length,
      n_peaks: peakRecords.length,
      n_tfs: tfRecords.length,
      n_edges: edgeRecords.length,
      n_snp_peaks: snpRecords.length,
      seed_genes: seedRecords,
      peaks: peakRecords,
      tfs: tfRecords,
      snps: snpRecords,
      edges: edgeRecords,
    };
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[c]));
  }
  function renderMessageMarkdown(s) {
    // Safe minimal markdown: escape first, then convert links + bold + newlines.
    let h = escapeHtml(s || "");
    h = h.replace(/\[([^\]]+)\]\((https?:\/\/[^)]+)\)/g,
      (_, t, u) => `<a href="${u}" target="_blank" rel="noopener">${t}</a>`);
    // Bare URLs → linkify (pubmed/doi/etc.)
    h = h.replace(/(?<![">])(https?:\/\/[^\s<>"']+)/g,
      (u) => `<a href="${u}" target="_blank" rel="noopener">${u}</a>`);
    // Convert graph-node mentions (rsids, peak IDs, gene/TF symbols) to
    // clickable spans → click focuses + highlights the node in cytoscape.
    h = inlineLinkifyNodes(h);
    h = h.replace(/\*\*([^*]+)\*\*/g, "<b>$1</b>");
    h = h.replace(/\n/g, "<br>");
    return h;
  }

  // Wrap any rsid / peak-ID / coordinate / gene / TF mention in chat prose
  // with `<span class="graph-node-link" data-id="...">…</span>` IF the
  // referenced node exists in the current Cytoscape graph. Re-uses the same
  // matching logic as extractMentionedNodeIds() for consistency.
  function inlineLinkifyNodes(html) {
    if (!S.cy) return html;
    const allNodes = new Set();
    S.cy.nodes().forEach((n) => allNodes.add(n.id()));
    const peakIds = [...allNodes].filter((id) => id.startsWith("chr"));

    function span(id, label, title) {
      const t = title ? ` title="${escapeHtml(title)}"` : "";
      return `<span class="graph-node-link" data-id="${escapeHtml(id)}"${t}>${label}</span>`;
    }

    // Don't linkify inside existing anchors — split by <a>...</a>.
    function transform(text) {
      // 1) Full peak IDs chrN:start-end
      text = text.replace(/\bchr[0-9XYMxym]+:\d+-\d+\b/g, (m) =>
        allNodes.has(m) ? span(m, m) : m);
      // 2) Single chrN:pos coords → resolve to containing peak
      text = text.replace(/\bchr[0-9XYMxym]+:\d+(?!-)/g, (m) => {
        if (allNodes.has(m)) return span(m, m);
        const p = m.match(/^(chr[0-9XYMxym]+):(\d+)/);
        if (!p) return m;
        const chr = p[1], pos = parseInt(p[2], 10);
        for (const pid of peakIds) {
          const pm = pid.match(/^(chr[0-9XYMxym]+):(\d+)-(\d+)$/);
          if (!pm || pm[1] !== chr) continue;
          const a = parseInt(pm[2], 10), b = parseInt(pm[3], 10);
          if (pos >= a && pos <= b) return span(pid, m, `peak ${pid}`);
        }
        return m;
      });
      // 3) rsids → look up peak via payload
      text = text.replace(/\brs\d+\b/g, (rs) => {
        for (const s of S.payload?.snps || []) {
          const otRs = s.annotations?.ot_variant?.rsIds || [];
          if ((s.rsid === rs || otRs.includes(rs)) && allNodes.has(s.peak)) {
            return span(s.peak, rs, `SNP ${rs} → peak ${s.peak}`);
          }
        }
        return rs;
      });
      // 4) Gene / TF symbols
      const ignore = new Set(["GWAS","SNP","DNA","RNA","AF","CAD","AVS","PMID",
        "DOI","URL","JSON","API","MHz","TBD","TFs","JACC","PMC","PNAS","NIH",
        "OK","ALL","MeSH","TF","FDR","HGNC","GTEx","VEP"]);
      text = text.replace(/\b[A-Z][A-Z0-9-]{1,}\b/g, (sym) => {
        if (ignore.has(sym)) return sym;
        return allNodes.has(sym) ? span(sym, sym) : sym;
      });
      return text;
    }

    // Split out any existing <a>...</a> sections so we don't double-wrap.
    const parts = html.split(/(<a\b[^>]*>.*?<\/a>)/g);
    return parts.map((p) => p.startsWith("<a") ? p : transform(p)).join("");
  }

  function focusGraphNode(id) {
    if (!S.cy) return;
    const n = S.cy.getElementById(id);
    if (!n.length) return;
    // Highlight in magenta + fit-on-node + open & pin the detail panel
    n.addClass("chat-highlight");
    setPinnedNode(n);             // → showDetailPanel(n) inside
    S.cy.animate({
      fit: { eles: n, padding: 180 },
      duration: 350, easing: "ease-out",
    });
  }
  function chatMsgEl(role, content, actions) {
    const div = document.createElement("div");
    div.className = `chat-msg ${role}`;
    // System / placeholder messages stay plain text; agent / user get markdown.
    if (role === "system") {
      div.textContent = content;
    } else {
      div.innerHTML = renderMessageMarkdown(content);
      // Click-to-focus on any node mentioned in the prose
      div.querySelectorAll(".graph-node-link").forEach((sp) => {
        sp.addEventListener("click", (e) => {
          e.preventDefault(); e.stopPropagation();
          focusGraphNode(sp.dataset.id);
        });
      });
    }
    if (actions && actions.length) {
      const ad = document.createElement("div");
      ad.className = "actions-applied";
      const labels = actions.map((a) => {
        const ids = a.args?.ids || (a.args?.from ? [a.args.from, a.args.to] : []);
        const idTxt = ids.length ? `: ${ids.slice(0, 6).join(", ")}${ids.length > 6 ? "…" : ""}` : "";
        return (a._auto ? "auto-highlight" : a.type) + idTxt;
      });
      ad.textContent = "▸ " + labels.join(" ｜ ");
      div.appendChild(ad);
    }
    return div;
  }

  function appendChatMessage(role, content, actions, citations, meta) {
    const body = $("chat-messages");
    body.appendChild(chatMsgEl(role, content, actions));
    if (Array.isArray(citations) && citations.length) {
      body.appendChild(citationsEl(citations, meta));
    } else if (meta?.grounding_count === 0 && role === "agent") {
      const note = document.createElement("div");
      note.className = "chat-msg system";
      note.textContent = "ℹ️ No PubMed citations grounded — agent prose is uncited inference.";
      body.appendChild(note);
    }
    body.scrollTop = body.scrollHeight;
  }

  function citationsEl(citations, meta) {
    const wrap = document.createElement("div");
    wrap.className = "chat-citations";
    const verified = citations.filter((c) => c.verified);
    const hallucinated = citations.filter((c) => c.hallucinated);
    const offGrnd = citations.filter((c) => c.off_grounding);
    let header = `<b>Citations (${verified.length} verified`;
    if (offGrnd.length) header += `, ${offGrnd.length} off-grounding`;
    if (hallucinated.length) header += `, ⚠️ ${hallucinated.length} hallucinated`;
    header += `)</b> <span class="meta">grounding=${meta?.grounding_count ?? "?"} from ${meta?.grounding_terms?.join(", ") || "—"}</span>`;
    wrap.innerHTML = `<div class="cit-header">${header}</div>`;

    const VISIBLE_DEFAULT = 5;
    function buildRow(c) {
      const div = document.createElement("div");
      div.className = "cit-row" + (c.verified ? "" : " cit-bad");
      const badge = c.hallucinated
        ? `<span class="cit-badge bad">⚠️ HALLUCINATED</span>`
        : c.off_grounding
        ? `<span class="cit-badge warn">off-grounding</span>`
        : `<span class="cit-badge ok">✓ verified</span>`;
      const linkTxt = c.pmid ? `PMID ${c.pmid}` : "(no PMID)";
      const url = c.url || (c.pmid ? `https://pubmed.ncbi.nlm.nih.gov/${c.pmid}/` : null);
      const link = url ? `<a href="${url}" target="_blank" rel="noopener">${linkTxt}</a>` : linkTxt;
      const meta2 = [c.year, c.journal, c.authors].filter(Boolean).join(" · ");
      div.innerHTML = `${badge} <b>${escapeHtml(c.title || "(no title)")}</b>`
        + `<div class="meta">${meta2} — ${link}</div>`
        + (c.key_finding ? `<div>${escapeHtml(c.key_finding)}</div>` : "");
      return div;
    }
    // Always show first VISIBLE_DEFAULT; hide overflow inside <details>
    const head = citations.slice(0, VISIBLE_DEFAULT);
    const tail = citations.slice(VISIBLE_DEFAULT);
    for (const c of head) wrap.appendChild(buildRow(c));
    if (tail.length) {
      const det = document.createElement("details");
      const sum = document.createElement("summary");
      sum.textContent = `Show ${tail.length} more`;
      sum.style.cssText = "cursor:pointer;font-size:11px;color:#1d4ed8;margin:4px 0";
      det.appendChild(sum);
      for (const c of tail) det.appendChild(buildRow(c));
      wrap.appendChild(det);
    }
    return wrap;
  }

  // Extract node IDs mentioned in agent prose. Returns IDs that EXIST in cy.
  // Handles: chr-region patterns (chrN:start-end or chrN:pos), rsids,
  // and uppercase gene/TF symbols that match an actual node.
  function extractMentionedNodeIds(text) {
    if (!text || !S.cy) return [];
    const all = new Set();
    S.cy.nodes().forEach((n) => all.add(n.id()));
    const peakIds = [...all].filter((id) => id.startsWith("chr"));
    const hits = new Set();

    // 1) Exact peak IDs (chrN:start-end)
    const peakRe = /\bchr[0-9XYMxym]+:\d+-\d+\b/g;
    for (const m of text.matchAll(peakRe)) {
      if (all.has(m[0])) hits.add(m[0]);
    }
    // 2) Single-coord mentions chrN:pos → find peak containing it
    const coordRe = /\bchr[0-9XYMxym]+:(\d+)(?!-)/g;
    for (const m of text.matchAll(coordRe)) {
      const fullMatch = m[0];
      if (all.has(fullMatch)) { hits.add(fullMatch); continue; }
      const chrPart = fullMatch.match(/^(chr[0-9XYMxym]+):(\d+)/);
      if (!chrPart) continue;
      const chr = chrPart[1];
      const pos = parseInt(chrPart[2], 10);
      for (const pid of peakIds) {
        const pm = pid.match(/^(chr[0-9XYMxym]+):(\d+)-(\d+)$/);
        if (!pm || pm[1] !== chr) continue;
        const start = parseInt(pm[2], 10), end = parseInt(pm[3], 10);
        if (pos >= start && pos <= end) hits.add(pid);
      }
    }
    // 3) rsids → look up payload to find their peak (rsid itself is not a node)
    const rsRe = /\brs\d+\b/g;
    for (const m of text.matchAll(rsRe)) {
      const rs = m[0];
      for (const s of S.payload?.snps || []) {
        const otRs = s.annotations?.ot_variant?.rsIds || [];
        if (s.rsid === rs || otRs.includes(rs)) {
          if (all.has(s.peak)) hits.add(s.peak);
        }
      }
    }
    // 4) Gene / TF symbols — 2+ uppercase chars, must exist in cy
    const symRe = /\b[A-Z][A-Z0-9-]{1,}\b/g;
    for (const m of text.matchAll(symRe)) {
      const sym = m[0];
      // Avoid common false positives
      if (["GWAS", "SNP", "DNA", "RNA", "AF", "CAD", "AVS", "PMID", "DOI",
           "URL", "JSON", "API", "MHz", "TBD", "TFs", "JACC",
           "PMC", "PNAS", "NIH", "OK"].includes(sym)) continue;
      if (all.has(sym)) hits.add(sym);
    }
    return [...hits];
  }

  async function sendChat(message) {
    if (S.chatPending) return;
    if (!message.trim()) return;
    S.chatPending = true;
    $("chat-send").disabled = true;
    appendChatMessage("user", message);
    S.chatHistory.push({ role: "user", content: message });

    // Inline "thinking" placeholder bubble in chat body with 1-sec ticker + bar
    const body = $("chat-messages");
    const ph = document.createElement("div");
    ph.className = "chat-msg agent chat-thinking";
    ph.innerHTML = `<div class="thinking-row">
        <span class="thinking-label">Agent thinking</span>
        <span class="thinking-elapsed">0s</span>
      </div>
      <div class="thinking-bar"><div class="thinking-bar-fill"></div></div>`;
    body.appendChild(ph);
    body.scrollTop = body.scrollHeight;
    const startTs = Date.now();
    const EXPECTED_SEC = 60;        // typical chat agent latency
    const tick = () => {
      const sec = Math.round((Date.now() - startTs) / 1000);
      const el = ph.querySelector(".thinking-elapsed");
      if (el) el.textContent = `${sec}s`;
      // Asymptotic fill: cap at 95% so the bar visibly settles before result arrives.
      const fill = Math.min(95, (sec / EXPECTED_SEC) * 90);
      const bar = ph.querySelector(".thinking-bar-fill");
      if (bar) bar.style.width = `${fill}%`;
    };
    tick();
    const timerId = setInterval(tick, 1000);
    const stopThinking = () => { clearInterval(timerId); ph.remove(); };

    // BYOK: load Anthropic key from localStorage; if missing, prompt the user
    const anthropicKey = localStorage.getItem("anthropic_api_key") || "";
    if (!anthropicKey || !anthropicKey.startsWith("sk-ant-")) {
      if (typeof window.openByokModal === "function") {
        stopThinking();
        appendChatMessage("system", "Chat needs your Anthropic API key. Opening setup…");
        window.openByokModal();
        return;
      } else {
        stopThinking();
        appendChatMessage("system", "Anthropic API key required (sk-ant-...). Set it via the Activate button.");
        return;
      }
    }

    try {
      const resp = await fetch(`${S.literatureBaseUrl}/chat`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "X-API-Key": anthropicKey,
        },
        body: JSON.stringify({
          message,
          graph_state: buildGraphState(),
          history: S.chatHistory.slice(-10),
          use_websearch: true,
        }),
      });
      if (resp.status === 401) {
        stopThinking();
        const data = await resp.json().catch(() => ({}));
        const why = data.detail || data.message || "Anthropic key was rejected.";
        appendChatMessage("system", `🔑 ${why} Re-enter your key.`);
        if (typeof window.openByokModal === "function") window.openByokModal();
        return;
      }
      if (resp.status === 429) {
        stopThinking();
        appendChatMessage("system", "⏳ Anthropic rate-limited your key. Wait a minute and try again.");
        return;
      }
      if (!resp.ok) {
        const txt = await resp.text();
        appendChatMessage("system", `Error ${resp.status}: ${txt.slice(0, 200)}`);
        return;
      }
      const data = await resp.json();
      const msg = data.message || (data._raw ? "(non-JSON response)" : "(empty)");
      const explicitActions = data.actions || [];

      // Pass 1: detect explicit disease / cell_state change requests
      let explicitDisease = null, explicitTargetCs = null;
      for (const a of explicitActions) {
        if (a.type === "set_disease") explicitDisease = a.args?.disease;
        if (a.type === "set_target_cs") explicitTargetCs = a.args?.cs;
      }
      // Auto-detect disease keywords from USER message (only) — clear intent
      const diseaseHints = {
        CAD: /\b(CAD|coronary artery disease|coronary disease|myocardial infarction)\b/i,
        AVS: /\b(AVS|aortic valve stenosis|aortic stenosis|calcific aortic)\b/i,
        AF:  /\b(AFib|atrial fibrillation)\b|^AF$|\bAF\b(?! disease)/i,
      };
      let autoDisease = null;
      for (const [d, re] of Object.entries(diseaseHints)) {
        if (re.test(message) && d !== S.diseaseId && S.manifest?.diseases[d]) {
          autoDisease = d; break;
        }
      }
      // Auto cell_state detection — STRICT mode.
      // Only fire when the USER message (not the agent reply) shows clear
      // switch intent. Naked mentions of a cs in passing should NOT trigger
      // a payload switch — that was causing the network to drift on every
      // comparison/contrast comment.
      const validTargets = S.manifest?.diseases[autoDisease || S.diseaseId]?.target_cs_list || [];
      const sortedCs = [...validTargets].sort((a, b) => b.length - a.length);
      // Informal cell_state hints (these are intent-strong by themselves)
      const informalCsHints = [
        [/\b(coronary\s+artery\s+SMC|coronary\s+SMC|vascular\s+SMC|SMC|smooth\s+muscle\s+cell)\b/i, "SmoothMuscleCellsCoronaryArtery"],
        [/\b(valve\s+fibroblast|VIC)\b/i, "ValveFibroblastsCardiacSkeleton"],
        [/\b(arterial\s+endothelial|EC\s+arterial)\b/i, "EndothelialCellsArterial"],
        [/\b(sleeve\s+cell|pulmonary\s+vein\s+sleeve|myocardial\s+sleeve)\b/i, "MyocardialSleeveCells"],
        [/\b(pacemaker\s+cell|SAN\s+cell)\b/i, "PacemakerCells"],
        [/\b(ventricular\s+cardiomyocyte|VCM|LV\s+myocyte)\b/i, "VentricularCardiomyocytesLeft"],
        [/\b(atrial\s+cardiomyocyte)\b/i, "AtrialCardiomyocytesLeft"],
        [/\b(activated\s+fibroblast|myofibroblast)\b/i, "FibroblastsActivated"],
      ];
      // Intent verb patterns. Each requires the cs name within ~3 words after.
      const csIntentTemplates = [
        // English: explicit "switch to", "show", "view", "load", "open" etc.
        (cs) => new RegExp(`\\b(switch to|move to|go to|load|open|view|look at|focus on|set target to|show( me)?|tell me about|about)\\W+(target\\s+)?${cs}\\b`, "i"),
        // English: tight "in/for/at/of" + cs (no slack — comparison phrases
        // like "compared to X" or "...in arrhythmogenesis where X..." won't match)
        (cs) => new RegExp(`\\b(in|for|at|of|on)\\s+${cs}\\b`, "i"),
        // cs at the start of a sentence / question
        (cs) => new RegExp(`(?:^|[.!?。\\n])\\s*${cs}\\b`, "i"),
        // Japanese particles that signal subject / topic / locative
        (cs) => new RegExp(`${cs}\\s*(について|の|を|では|で|に)`, "i"),
        // Japanese imperative / question pattern targeting the cs
        (cs) => new RegExp(`(見せて|表示して|教えて|に切替|に変更|に移動|を選んで|を開いて)[^。\\n]{0,15}${cs}`, "i"),
      ];
      let autoTarget = null;
      let autoTargetReason = null;
      // 1) Strict cs detection in USER message only
      for (const cs of sortedCs) {
        if (cs === S.targetCs) continue;
        for (const tmpl of csIntentTemplates) {
          if (tmpl(cs).test(message)) {
            autoTarget = cs;
            autoTargetReason = `user intent verb + cs name`;
            break;
          }
        }
        if (autoTarget) break;
      }
      // 2) Informal keyword hints (intent-strong by themselves)
      if (!autoTarget) {
        for (const [re, preferred] of informalCsHints) {
          if (!re.test(message)) continue;
          if (validTargets.includes(preferred) && preferred !== S.targetCs) {
            autoTarget = preferred;
            autoTargetReason = `informal hint match`;
            break;
          }
          const family = preferred.replace(/Cells.*$|Left.*$|Right.*$|Coronary.*$|Capillary.*$|Cardiac.*$|Pos$|Neg$/, "");
          const match = sortedCs.find((c) => c.startsWith(family) && c !== S.targetCs);
          if (match) { autoTarget = match; autoTargetReason = `informal hint family fallback`; break; }
        }
      }
      if (autoTarget) console.log("[chat auto-cs]", autoTarget, "(", autoTargetReason, ")");
      // Build initial action list (disease + target/ref switches go first)
      const allActions = [...explicitActions];
      if (autoDisease && !explicitDisease) {
        allActions.push({ type: "set_disease", args: { disease: autoDisease }, _auto: true });
      }
      if (autoTarget && !explicitTargetCs) {
        allActions.push({ type: "set_target_cs", args: { cs: autoTarget }, _auto: true });
      }
      // Run disease/target/ref switch actions FIRST and wait — subsequent
      // highlights will be evaluated against the newly-rebuilt subgraph.
      const switchActions = allActions.filter((a) =>
        a.type === "set_disease" || a.type === "set_target_cs" || a.type === "set_ref_cs"
      ).sort((x, y) => ({ set_disease: 0, set_target_cs: 1, set_ref_cs: 2 }[x.type] -
                          { set_disease: 0, set_target_cs: 1, set_ref_cs: 2 }[y.type]));
      for (const a of switchActions) {
        const fn = ACTION_HANDLERS[a.type];
        try { await fn(a.args || {}); } catch (e) { console.error(e); }
      }
      // NOW extract IDs (against fresh subgraph) and add to remaining actions
      const explicitIds = new Set();
      for (const a of explicitActions) {
        if (a.type === "highlight_nodes") for (const id of (a.args?.ids || [])) explicitIds.add(id);
        if (a.type === "highlight_path") { if (a.args?.from) explicitIds.add(a.args.from); if (a.args?.to) explicitIds.add(a.args.to); }
      }
      const autoIds = extractMentionedNodeIds(msg).filter((id) => !explicitIds.has(id));
      const remainingActions = allActions.filter((a) =>
        a.type !== "set_disease" && a.type !== "set_target_cs" && a.type !== "set_ref_cs");
      if (autoIds.length) {
        remainingActions.push({ type: "highlight_nodes", args: { ids: autoIds }, _auto: true });
      }
      appendChatMessage("agent", msg, allActions, data.citations || [],
                        { grounding_count: data.grounding_count,
                          grounding_terms: data.grounding_terms });
      S.chatHistory.push({ role: "agent", content: msg, actions: allActions });
      await applyChatActions(remainingActions);
      if (data._raw) {
        const raw = document.createElement("div");
        raw.className = "raw";
        raw.textContent = "raw: " + data._raw.slice(0, 800);
        $("chat-messages").lastElementChild.appendChild(raw);
      }
    } catch (err) {
      appendChatMessage("system", "Backend unreachable — start with: python scripts/snp2cell_hub_backend.py");
    } finally {
      stopThinking();
      S.chatPending = false;
      $("chat-send").disabled = false;
    }
  }

  // ---------- action dispatcher ----------
  const ACTION_HANDLERS = {
    highlight_nodes: (args) => {
      if (!S.cy) return [];
      const ids = (args && args.ids) || [];
      const ok = [];
      S.cy.batch(() => {
        ids.forEach((id) => {
          const n = S.cy.getElementById(id);
          if (n.length) { n.addClass("chat-highlight"); ok.push(id); }
        });
      });
      return ok;
    },
    highlight_path: (args) => {
      if (!S.cy) return [];
      const from = S.cy.getElementById(args.from);
      const to = S.cy.getElementById(args.to);
      if (!from.length || !to.length) return [];
      const r = S.cy.elements().aStar({ root: from, goal: to, directed: false });
      if (!r.found) return [];
      r.path.addClass("chat-highlight");
      return r.path.nodes().map((n) => n.id());
    },
    highlight_cell_type_edges: (args) => {
      if (!S.cy) return [];
      const ct = args && args.cell_type;
      if (!ct) return [];
      const matched = [];
      S.cy.batch(() => {
        S.cy.edges().forEach((e) => {
          const cts = (e.data("cts") || "").split(", ");
          if (cts.includes(ct)) { e.addClass("chat-highlight"); matched.push(e.id()); }
        });
      });
      return matched;
    },
    set_view: (args) => {
      const a = args || {};
      const out = [];
      if (Number.isFinite(a.top_n)) {
        S.topN = Math.max(10, Math.min(200, a.top_n));
        $("topn-slider").value = S.topN; $("topn-label").textContent = S.topN;
        out.push(`top_n=${S.topN}`);
      }
      if (Number.isFinite(a.peaks_per_seed)) {
        S.peaksPerSeed = Math.max(1, Math.min(10, a.peaks_per_seed));
        $("peaks-slider").value = S.peaksPerSeed; $("peaks-label").textContent = S.peaksPerSeed;
        out.push(`peaks_per_seed=${S.peaksPerSeed}`);
      }
      if (Number.isFinite(a.tfs_per_peak)) {
        S.tfsPerPeak = Math.max(0, Math.min(10, a.tfs_per_peak));
        $("tfs-slider").value = S.tfsPerPeak; $("tfs-label").textContent = S.tfsPerPeak;
        out.push(`tfs_per_peak=${S.tfsPerPeak}`);
      }
      if (Array.isArray(a.edge_highlight_cell_types)) {
        S.edgeHighlightCts = a.edge_highlight_cell_types;
        // sync select element
        for (const opt of $("edge-ct-select").options) opt.selected = S.edgeHighlightCts.includes(opt.value);
        out.push(`edge_highlight=[${S.edgeHighlightCts.join(",")}]`);
      }
      rebuild();
      return out;
    },
    focus_on_seed: (args) => {
      if (!S.cy) return [];
      const id = args && args.id;
      if (!id) return [];
      const n = S.cy.getElementById(id);
      if (!n.length) return [];
      const keep = new Set([id, ...n.neighborhood("node").map((x) => x.id())]);
      S.cy.batch(() => {
        S.cy.nodes().forEach((nd) => {
          if (!keep.has(nd.id())) nd.addClass("faded");
          else nd.addClass("chat-highlight");
        });
        S.cy.edges().forEach((e) => {
          if (!keep.has(e.source().id()) || !keep.has(e.target().id())) e.addClass("faded");
        });
      });
      S.cy.fit(S.cy.elements(".chat-highlight"), 60);
      return [id];
    },
    search_pubmed: (args) => {
      const a = args || {};
      if (!a.type || !a.key) return [];
      searchLiterature(a.type, a.key, a.mode || "agent");
      return [`${a.type}/${a.key}`];
    },
    set_disease: async (args) => {
      const d = args?.disease;
      if (!d) return [];
      if (!S.manifest?.diseases[d]) {
        console.warn(`set_disease: unknown disease ${d}`);
        return [];
      }
      if (d === S.diseaseId) return [`already on ${d}`];
      await selectDisease(d);
      return [`disease=${d}`];
    },
    set_target_cs: async (args) => {
      const cs = args?.cs;
      if (!cs) return [];
      const valid = S.manifest?.diseases[S.diseaseId]?.target_cs_list || [];
      if (!valid.includes(cs)) {
        console.warn(`set_target_cs: ${cs} not in built target_cs_list (need rebuild)`);
        return [];
      }
      if (cs === S.targetCs) return [`already on ${cs}`];
      await selectTargetCs(cs);
      return [`target_cs=${cs}`];
    },
    set_ref_cs: (args) => {
      const list = (args?.cs_list || []).filter((c) => c !== S.targetCs
                    && (S.deZ?.csIdx?.has(c) ?? false));
      if (list.length < 2) {
        console.warn("set_ref_cs needs ≥2 valid cs (got", list.length, ")");
        return [];
      }
      S.refMode = "selected";
      S.refCsList = list;
      $("ref-mode").value = "selected";
      $("ref-picker").style.display = "block";
      populateRefPicker();
      debouncedRecomputeSpec();
      return [`ref_cs=[${list.join(",")}]`];
    },
    reset_highlights: () => {
      if (!S.cy) return [];
      S.cy.elements().removeClass("chat-highlight faded");
      return ["all"];
    },
  };

  async function applyChatActions(actions) {
    if (!actions || !actions.length) return;
    // set_target_cs first (it reloads the payload, which would invalidate
    // any subsequent highlight_nodes IDs from the OLD subgraph)
    const order = (a) => ({
      set_disease: 0, set_target_cs: 1, set_ref_cs: 2, set_view: 3,
      reset_highlights: 4, focus_on_seed: 5, highlight_path: 6,
      highlight_nodes: 7, highlight_cell_type_edges: 8, search_pubmed: 9,
    }[a.type] ?? 99);
    const sorted = [...actions].sort((x, y) => order(x) - order(y));
    for (const a of sorted) {
      const fn = ACTION_HANDLERS[a.type];
      if (!fn) {
        console.warn("Unknown chat action:", a.type);
        continue;
      }
      try {
        const result = await fn(a.args || {});
        console.log("[chat action]", a.type, a.args, "→", result);
      } catch (err) {
        console.error("Action error:", err);
      }
    }
  }

  function openChat() {
    const p = $("chat-panel");
    p.style.display = "flex";
    p.classList.remove("minimized");
    $("chat-input").focus();
  }
  function openChatAndPrefill(kind, key) {
    openChat();
    const tpl = {
      snp:  `Tell me about ${key}: literature in the context of ${S.diseaseId}, links to genes/TFs in this network, and highlight on the graph.`,
      gene: `What is the role of ${key} in ${S.targetCs} / ${S.diseaseId}? Cite relevant papers and highlight ${key} + immediate neighbors.`,
      tf:   `Summarise the role of TF ${key} in ${S.targetCs}: downstream targets in the current graph, supporting literature, and highlight TBX-style path.`,
      peak: `Tell me about peak ${key}: associated SNPs, regulated genes in this network, and relevant papers.`,
    };
    const text = tpl[kind] || `Tell me about ${key}.`;
    $("chat-input").value = text;
    // Auto-resize and select for easy editing
    $("chat-input").setSelectionRange(text.length, text.length);
  }

  function setChatMinimized(min) {
    const p = $("chat-panel");
    p.classList.toggle("minimized", min);
    $("chat-close-btn").textContent = min ? "▢" : "−";
    $("chat-close-btn").title = min ? "Expand" : "Minimize";
    if (!min) $("chat-input").focus();
  }
  function wireChatPanel() {
    $("chat-toggle-btn").addEventListener("click", () => {
      // Cycle: hidden → expanded → minimized → expanded → ...
      const p = $("chat-panel");
      if (p.style.display === "none") {
        p.style.display = "flex"; setChatMinimized(false);
      } else if (p.classList.contains("minimized")) {
        setChatMinimized(false);
      } else {
        setChatMinimized(true);
      }
    });
    // Minimize button on header
    $("chat-close-btn").addEventListener("click", (e) => {
      e.stopPropagation();
      const p = $("chat-panel");
      setChatMinimized(!p.classList.contains("minimized"));
    });
    // Click on minimized header → restore
    $("chat-panel").querySelector(".chat-header").addEventListener("click", (e) => {
      if (e.target.id === "chat-close-btn" || e.target.id === "chat-clear") return;
      if ($("chat-panel").classList.contains("minimized")) setChatMinimized(false);
    });
    $("chat-clear").addEventListener("click", () => {
      S.chatHistory = [];
      $("chat-messages").innerHTML = "";
      appendChatMessage("system", "History cleared.");
    });
    $("chat-send").addEventListener("click", () => {
      const v = $("chat-input").value;
      $("chat-input").value = "";
      sendChat(v);
    });
    $("chat-input").addEventListener("keydown", (e) => {
      if (e.key === "Enter" && (e.metaKey || e.ctrlKey)) {
        e.preventDefault();
        $("chat-send").click();
      }
    });
    // Drag header
    const panel = $("chat-panel");
    const header = panel.querySelector(".chat-header");
    let drag = null;
    header.addEventListener("mousedown", (e) => {
      if (e.target.tagName === "BUTTON" || e.target.classList.contains("chat-close")) return;
      const r = panel.getBoundingClientRect();
      drag = { dx: e.clientX - r.left, dy: e.clientY - r.top };
    });
    document.addEventListener("mousemove", (e) => {
      if (!drag) return;
      panel.style.left = `${e.clientX - drag.dx}px`;
      panel.style.top = `${e.clientY - drag.dy}px`;
      panel.style.right = "auto"; panel.style.bottom = "auto";
    });
    document.addEventListener("mouseup", () => { drag = null; });
  }

  // ---------- AG1 (AlphaGenome) variant-effect panel ----------
  let _ag1Panel = null;
  let _ag1Polls = new Map();
  const HUB_BASE = "http://127.0.0.1:8765/hub";   // for static cache files

  function ensureAg1Panel() {
    if (_ag1Panel) return _ag1Panel;
    const p = document.createElement("div");
    p.className = "ag1-panel";
    p.style.display = "none";
    p.innerHTML = `<div class="ag1-header">
        <span>🧬 AG1 variant effect <span id="ag1-title" class="meta"></span></span>
        <span class="ag1-actions">
          <button id="ag1-reset-pos" title="Reset position">⤢</button>
          <span class="chat-close" id="ag1-close" title="Close (Esc)">×</span>
        </span>
      </div>
      <div class="ag1-body" id="ag1-body"></div>`;
    document.body.appendChild(p);
    p.querySelector("#ag1-close").addEventListener("click", () => p.style.display = "none");
    p.querySelector("#ag1-reset-pos").addEventListener("click", () => resetAg1Position(p));
    // drag header
    const header = p.querySelector(".ag1-header");
    let drag = null;
    header.addEventListener("mousedown", (e) => {
      if (e.target.id === "ag1-close") return;
      const r = p.getBoundingClientRect();
      drag = { dx: e.clientX - r.left, dy: e.clientY - r.top };
    });
    document.addEventListener("mousemove", (e) => {
      if (!drag) return;
      // Clamp so header (with close button) stays in viewport
      const rect = p.getBoundingClientRect();
      const minTop = 0, minLeft = 24 - rect.width;     // leave at least 24px visible
      const maxTop = window.innerHeight - 32;          // keep header visible
      const maxLeft = window.innerWidth - 24;
      let nx = e.clientX - drag.dx;
      let ny = e.clientY - drag.dy;
      nx = Math.max(minLeft, Math.min(nx, maxLeft));
      ny = Math.max(minTop, Math.min(ny, maxTop));
      p.style.left = `${nx}px`;
      p.style.top = `${ny}px`;
      p.style.right = "auto"; p.style.bottom = "auto";
    });
    document.addEventListener("mouseup", () => { drag = null; });
    // ESC key closes the panel
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape" && p.style.display === "block") {
        p.style.display = "none";
      }
    });
    _ag1Panel = p;
    return p;
  }

  function resetAg1Position(p) {
    p.style.top = "60px"; p.style.right = "60px";
    p.style.left = "auto"; p.style.bottom = "auto";
  }

  function populateAg1CsSelect(filterStr) {
    const sel = $("ag1-cs-input");
    if (!sel) return;
    const groups = S.lineage?.groups || {};
    const nCells = S.lineage?.n_cells || {};
    const q = (filterStr || "").trim().toLowerCase();
    // Preserve current selections across re-filter
    const prevSelected = new Set(Array.from(sel.selectedOptions || []).map((o) => o.value));
    if (!prevSelected.size && S.targetCs) prevSelected.add(S.targetCs);

    sel.innerHTML = "";
    const types = Object.keys(groups).sort();
    let total = 0;
    for (const t of types) {
      const inGroup = groups[t].filter((cs) => !q || cs.toLowerCase().includes(q));
      if (!inGroup.length) continue;
      const og = document.createElement("optgroup");
      og.label = t;
      for (const cs of inGroup.sort()) {
        const opt = document.createElement("option");
        opt.value = cs;
        const n = nCells[cs];
        opt.textContent = `${cs}${n != null ? ` (n=${n.toLocaleString()})` : ""}`;
        if (prevSelected.has(cs)) opt.selected = true;
        og.appendChild(opt);
        total++;
      }
      sel.appendChild(og);
    }
    const countEl = $("ag1-cs-count");
    if (countEl) {
      const nSel = Array.from(sel.selectedOptions || []).length;
      countEl.textContent = `${nSel} selected · ${total}${q ? ` / ${Object.values(groups).flat().length} match` : ""}`;
    }
  }

  // Wire after panel is populated
  function _ag1WireSelectHandlers() {
    const filter = $("ag1-cs-filter");
    const sel = $("ag1-cs-input");
    if (!filter || !sel) return;
    filter.addEventListener("input", (e) => populateAg1CsSelect(e.target.value));
    sel.addEventListener("change", () => populateAg1CsSelect($("ag1-cs-filter").value));
  }

  async function openAg1Panel(variantId, rsid) {
    const p = ensureAg1Panel();
    p.style.display = "block";
    // If previously dragged off-screen, snap back into view on open
    const rect = p.getBoundingClientRect();
    if (rect.top < 0 || rect.left < 0 || rect.left > window.innerWidth - 40
        || rect.top > window.innerHeight - 40) {
      resetAg1Position(p);
    }
    $("ag1-title")?.replaceChildren?.(document.createTextNode(`${rsid} · ${variantId}`));
    const body = $("ag1-body");
    // Probe current cache state
    let status = null;
    try {
      const r = await fetch(`${S.literatureBaseUrl}/ag1/status/${variantId}`);
      if (r.ok) status = await r.json();
    } catch (_) {}
    if (status?.status === "done") {
      renderAg1Result(variantId, rsid, status);
      return;
    }
    body.innerHTML = `<div class="meta">Variant <code>${variantId}</code> not yet scored.</div>
      <div class="ag1-warn">⚠ AG1 inference is heavy — about <b>6-8 minutes on CPU</b>. Result is cached.</div>
      <button class="ag1-start-btn" data-variant="${variantId}">🧬 Start scoring</button>`;
    body.querySelector(".ag1-start-btn").addEventListener("click", async () => {
      body.innerHTML = `<div class="ag1-progress">
          <div class="meta">Spawning AG1 subprocess…</div>
          <div class="ag1-elapsed" id="ag1-elapsed-${variantId}">0s · est. ~480s</div>
          <div class="thinking-bar"><div class="thinking-bar-fill" id="ag1-bar-${variantId}"></div></div>
          <div class="meta" style="margin-top:6px">Will auto-refresh when done. You can close this panel; the run continues in background.</div>
        </div>`;
      try {
        await fetch(`${S.literatureBaseUrl}/ag1/score/${variantId}`, {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({}),
        });
      } catch (e) {
        body.innerHTML = `<div class="ag1-warn">Backend unreachable.</div>`;
        return;
      }
      pollAg1(variantId, rsid);
    });
  }

  function pollAg1(variantId, rsid) {
    if (_ag1Polls.has(variantId)) return;
    const startTs = Date.now();
    const EXPECTED_SEC = 480;
    // 1-second tick — pure UI (no network)
    const tick = () => {
      const sec = Math.round((Date.now() - startTs) / 1000);
      const el = document.getElementById(`ag1-elapsed-${variantId}`);
      if (el) el.textContent = `${sec}s · est. ~${EXPECTED_SEC}s`;
      const bar = document.getElementById(`ag1-bar-${variantId}`);
      if (bar) bar.style.width = `${Math.min(95, (sec / EXPECTED_SEC) * 90)}%`;
    };
    tick();
    const uiTick = setInterval(tick, 1000);
    // 8-second poll — actual backend status fetch
    const t = setInterval(async () => {
      try {
        const r = await fetch(`${S.literatureBaseUrl}/ag1/status/${variantId}`);
        if (!r.ok) return;
        const s = await r.json();
        if (s.status === "done") {
          clearInterval(t); clearInterval(uiTick); _ag1Polls.delete(variantId);
          if (_ag1Panel && _ag1Panel.style.display === "block"
              && $("ag1-title").textContent.includes(variantId)) {
            renderAg1Result(variantId, rsid, s);
          }
        } else if (s.status === "error") {
          clearInterval(t); clearInterval(uiTick); _ag1Polls.delete(variantId);
          if (_ag1Panel) {
            $("ag1-body").innerHTML = `<div class="ag1-warn">AG1 error: ${s.error || "(no detail)"}</div>`;
          }
        }
      } catch (_) { /* keep polling */ }
    }, 8000);
    _ag1Polls.set(variantId, t);
  }

  function renderAg1Result(variantId, rsid, status) {
    const body = $("ag1-body");
    const atac = status.atac?.top || [];
    const rna  = status.rna?.top  || [];
    const fmtScore = (v) => (v == null ? "—" : v.toFixed(3));
    function renderTable(rows, title, kind) {
      if (!rows.length) return `<div class="meta">No ${title} hits.</div>`;
      return `<div class="ag1-section"><b>${title}</b> — ${rows.length} rows (sorted |effect| desc, scroll for more)</div>
        <div class="ag1-table-scroll">
          <table class="ag1-table">
            <thead><tr><th>cell_state</th><th>${kind === "rna" ? "gene" : ""}</th><th>raw</th><th>|raw|</th></tr></thead>
            <tbody>
            ${rows.map((r) => `<tr>
                <td>${r.cell_state || ""}</td>
                <td>${kind === "rna" ? (r.gene_name || "") : ""}</td>
                <td class="${r.raw_score >= 0 ? "pos" : "neg"}">${fmtScore(r.raw_score)}</td>
                <td>${fmtScore(r.abs_score)}</td>
              </tr>`).join("")}
            </tbody>
          </table>
        </div>`;
    }
    const csvLink = (rel) => rel ? `<a target="_blank" href="${HUB_BASE}/${rel}">${rel.split('/').pop()}</a>` : "";
    body.innerHTML = `
      <div class="meta">Variant <code>${variantId}</code> — cached. ATAC ${status.atac?.n_rows || 0} rows · RNA ${status.rna?.n_rows || 0} rows.</div>
      <div class="meta">CSV: ${csvLink(status.atac_csv)} · ${csvLink(status.rna_csv)}</div>
      ${renderTable(atac, "ATAC", "atac")}
      ${renderTable(rna,  "RNA",  "rna")}
      <div class="ag1-section"><b>REF / ALT / Δ track</b></div>
      <div class="control">
        <label>Cell states <span class="meta">(Cmd/Ctrl+click for multiple)</span></label>
        <input id="ag1-cs-filter" type="text" placeholder="filter…" style="margin-bottom:4px">
        <select id="ag1-cs-input" multiple size="9" style="height:auto;width:100%"></select>
        <div class="meta" id="ag1-cs-count" style="margin-top:2px"></div>
      </div>
      <div class="control">
        <label>Window (bp)</label>
        <input id="ag1-window-input" type="number" value="400000" step="50000" min="50000" max="1500000">
      </div>
      <button class="ag1-render-btn" data-variant="${variantId}">Render track →</button>
      <div id="ag1-render-out"></div>
    `;
    // Populate cell-state multi-select grouped by cell_type + wire filter
    populateAg1CsSelect();
    _ag1WireSelectHandlers();
    body.querySelector(".ag1-render-btn").addEventListener("click", async () => {
      const sel = $("ag1-cs-input");
      const cs = Array.from(sel.selectedOptions).map((o) => o.value);
      if (!cs.length) {
        $("ag1-render-out").innerHTML = `<div class="ag1-warn">Pick at least one cell state.</div>`;
        return;
      }
      const win = parseInt($("ag1-window-input").value || 400000);
      const out = $("ag1-render-out");
      out.innerHTML = `<div class="meta">Rendering (will take ~10-30 s)…</div>`;
      try {
        const r = await fetch(`${S.literatureBaseUrl}/ag1/render/${variantId}`, {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ cell_states: cs, window: win }),
        });
        const data = await r.json();
        if (data.error) {
          out.innerHTML = `<div class="ag1-warn">Render error: ${data.error}</div>`;
          return;
        }
        out.innerHTML = `<div><a target="_blank" href="${HUB_BASE}/${data.png}">${data.png.split('/').pop()}</a></div>
          <img src="${HUB_BASE}/${data.png}?t=${Date.now()}" class="ag1-img" alt="track">`;
      } catch (e) {
        out.innerHTML = `<div class="ag1-warn">Backend unreachable.</div>`;
      }
    });
  }

  // ---------- bootstrap ----------
  document.addEventListener("DOMContentLoaded", () => { init().then(wireChatPanel); });
})();
