// snp2cell hub: client-side specificity recompute.
// Mirrors snp2cell_figures.specificity.compute_specificity_padj:
//   z = (target - mean(ref)) / std(ref, ddof=1)
//   p = 1 - tcdf(z, df=n_ref - 1)
//   padj = BH(p) across all genes
//
// All math here is plain JS — no dependencies. Used by viewer.js on every
// ref_cs toggle.

(function () {
  "use strict";

  // ---- Gamma / Beta (Numerical Recipes lite) ----
  function gammln(x) {
    const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
                 -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x;
    let tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) {
      y += 1;
      ser += cof[j] / y;
    }
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  function betacf(a, b, x) {
    const FPMIN = 1e-30;
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1;
    let d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 200; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d;
      if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c;
      if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d;
      if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c;
      if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 3e-9) break;
    }
    return h;
  }

  function betai(a, b, x) {
    if (Number.isNaN(x) || x < 0 || x > 1) return NaN;
    if (x === 0 || x === 1) return x;
    const bt = Math.exp(gammln(a + b) - gammln(a) - gammln(b)
                        + a * Math.log(x) + b * Math.log(1 - x));
    if (x < (a + 1) / (a + b + 2)) return bt * betacf(a, b, x) / a;
    return 1 - bt * betacf(b, a, 1 - x) / b;
  }

  // Student's t cumulative distribution P(T <= t | df)
  function tCdf(t, df) {
    if (df <= 0 || Number.isNaN(t)) return NaN;
    if (!Number.isFinite(t)) return t > 0 ? 1 : 0;
    const x = df / (df + t * t);
    const ix = betai(df / 2, 0.5, x);
    return t >= 0 ? 1 - 0.5 * ix : 0.5 * ix;
  }

  // Upper-tail p value for one-sided greater test
  function tSf(t, df) {
    if (df <= 0 || Number.isNaN(t)) return NaN;
    if (!Number.isFinite(t)) return t > 0 ? 0 : 1;
    return 1 - tCdf(t, df);
  }

  // Benjamini-Hochberg adjustment, equivalent to statsmodels multipletests(method="fdr_bh")
  function bhAdjust(pvals) {
    const n = pvals.length;
    const idx = new Int32Array(n);
    for (let i = 0; i < n; i++) idx[i] = i;
    // Sort indices by p ascending; NaN treated as 1
    idx.sort((a, b) => {
      const pa = Number.isFinite(pvals[a]) ? pvals[a] : 1;
      const pb = Number.isFinite(pvals[b]) ? pvals[b] : 1;
      return pa - pb;
    });
    const adj = new Float64Array(n);
    let prev = 1;
    for (let k = n - 1; k >= 0; k--) {
      const rank = k + 1;
      const p = Number.isFinite(pvals[idx[k]]) ? pvals[idx[k]] : 1;
      const v = Math.min(prev, p * n / rank);
      prev = v;
      adj[idx[k]] = Math.max(0, Math.min(1, v));
    }
    return adj;
  }

  // Core call: returns {z, p, padj, passesFdr} all Float64Array(n_genes).
  // Inputs:
  //   matrix       : Float32Array of shape n_genes × n_cs (row-major)
  //   nGenes, nCs  : matrix dims
  //   targetIdx    : column index for target_cs
  //   refIdx       : Int32Array of column indices for reference cs (≥2)
  //   fdrThresh    : threshold for passesFdr mask
  function computeSpecificity(matrix, nGenes, nCs, targetIdx, refIdx, fdrThresh) {
    if (refIdx.length < 2) {
      throw new Error(`reference_cs_list needs ≥2 entries (got ${refIdx.length})`);
    }
    const df = refIdx.length - 1;
    const pvals = new Float64Array(nGenes);
    const zs = new Float64Array(nGenes);

    for (let g = 0; g < nGenes; g++) {
      const base = g * nCs;
      const target = matrix[base + targetIdx];
      // Welford-ish for mean & std (ddof=1)
      let sum = 0;
      for (let k = 0; k < refIdx.length; k++) sum += matrix[base + refIdx[k]];
      const mean = sum / refIdx.length;
      let ss = 0;
      for (let k = 0; k < refIdx.length; k++) {
        const v = matrix[base + refIdx[k]] - mean;
        ss += v * v;
      }
      const sd = Math.sqrt(ss / df);
      if (!Number.isFinite(sd) || sd === 0) {
        zs[g] = NaN; pvals[g] = NaN;
        continue;
      }
      const z = (target - mean) / sd;
      zs[g] = z;
      pvals[g] = tSf(z, df);   // upper tail (greater)
    }

    const padj = bhAdjust(pvals);
    const passesFdr = new Uint8Array(nGenes);
    for (let g = 0; g < nGenes; g++) {
      passesFdr[g] = (Number.isFinite(padj[g]) && padj[g] < fdrThresh) ? 1 : 0;
    }
    return { z: zs, p: pvals, padj, passesFdr, df };
  }

  window.SpecificityMath = { tCdf, tSf, bhAdjust, computeSpecificity, gammln };
})();
