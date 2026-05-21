# CardioNav — User Guide

*Last updated: 2026-05-21*

---

## 1. What is CardioNav?

CardioNav is an interactive viewer for **pre-computed cardiac gene
regulatory networks (eGRNs)** that integrate GWAS credible-set fine-mapping,
single-cell RNA-seq, single-cell ATAC-seq, and Open Targets annotations.

Coverage today:

- **3 diseases**: atrial fibrillation (AF), coronary artery disease (CAD),
  aortic valve stenosis (AVS)
- **~50 cell types / cell states per disease** — atrial cardiomyocytes,
  endothelial cells (arterial / venous / capillary subsets), fibroblasts,
  valve interstitial cells, immune cell subsets, conduction system, etc.
- **PubMed-grounded chat agent** that can answer questions about the
  current view, switch the viewer for you, and cite verified PMIDs.

The networks were built from a `snp2cell` + `SCENIC+` + Open Targets v26.03
pipeline. **CardioNav is a viewer, not a pipeline runner** — the underlying
eGRNs are snapshots and cannot be regenerated from inside the app.

### What CardioNav is good for

- Exploring **regulatory hypotheses** at GWAS loci (which TFs likely drive
  which seed genes in which cell state).
- Finding **cell-state specific signal** — which loci are active where.
- **Cross-locus pivoting** — does my SNP of interest appear in any other
  disease / cell-state combo?
- **Literature-anchored brainstorming** — quick LLM-mediated lookup with
  inline PubMed citations.

### What CardioNav is not

- Not a clinical-decision tool.
- Not an on-demand variant-effect predictor — the AlphaGenome panel covers
  only 15 representative variants (beta).
- Not a primary source of truth for any single biological claim — every
  network edge is a **predicted** regulatory link, and every citation
  returned by the chat agent should be read by a human before quoting.

---

## 2. Quickstart (5 minutes)

### 2.1  Log in

The site lives at the URL your collaborator sent you. Your browser will
prompt for a username + password (HTTP Basic Auth). The credentials are
shared by the PI — there is no per-user account in this preview.

### 2.2  Pick a disease and a cell state

Top bar, left side:

- **Disease dropdown** → AF / CAD / AVS
- **Target cell-state dropdown** → e.g. `AtrialCardiomyocytesLeft`
- (Optional) **Reference cell-state** → controls the DE-z column used to
  colour seed genes

The graph rebuilds automatically. Loading the **DE-z matrix** the first
time takes a few seconds (it's a ~5 MB gzipped matrix); subsequent visits
are cache-served.

### 2.3  Read the graph

Concentric layout, three rings:

| Ring | What it represents |
|---|---|
| **Outer** | Seed genes — disease-prioritised by GWAS-z |
| **Middle** | ATAC peaks — regulatory regions overlapping credible-set SNPs |
| **Inner** | Transcription factors — predicted to bind those peaks |

Edges:
- **TF → peak** (inner → middle)
- **peak → gene** (middle → outer)

Edges are coloured by the supporting cell-type set. Use the top-N sliders
(top right) to thin the network for clarity.

### 2.4  Inspect a SNP

Hover over any peak or click a SNP node → a tooltip card opens with:

- **rsid**, chr/pos/ref/alt
- **PIP** (posterior inclusion probability from Open Targets fine-mapping)
- **gnomAD allele frequency**
- **trait** (study label) and **`study_id`** (specific OT credibleSet)
- **VEP most-severe consequence**
- **L2G top gene** (probabilistic SNP → gene mapping)
- Link to PubMed / Open Targets variant page
- 🧬 button to open the AlphaGenome panel (beta, only some variants)

### 2.5  Set up the chat agent

Click 🔑 (top right). Paste either:

- An **Anthropic key** (`sk-ant-…`) — better instruction-following, but
  free tier is non-existent; you need a paid Anthropic plan
- A **Google Gemini key** (`AIza…`) — generous free tier
  (`gemini-2.0-flash`: 15 RPM, 1M tokens/min, 1500 RPD)

The key is stored only in your browser's localStorage. It is sent per
request as an `X-API-Key` header over HTTPS and never persisted on the
server. See [How API keys are protected](#api-key-protection) below.

Ask a question in the 💬 chat panel (right side). The agent can:

- Answer questions about the current view
- Switch the viewer ("show me TBX5 in atrial cardiomyocytes")
- Cite PubMed PMIDs inline — click a citation card to read

---

## 3. Core concepts primer

Five terms that experienced GWAS folks sometimes still misinterpret.

### 3.1  eGRN — enhancer Gene Regulatory Network

A **disease- and cell-state-specific** graph of *predicted* regulatory
relationships: which TFs are likely to bind which open chromatin regions
(ATAC peaks), and which genes those peaks are likely to regulate.

> **eGRN edges are hypotheses derived from sequence + chromatin + expression
> data — not experimentally validated interactions.**

A high-confidence edge in CardioNav means the pipeline saw consistent
evidence (TF motif in an accessible peak that loops to a gene whose
expression correlates with the TF's). It does NOT mean a ChIP-seq or
CRISPRi experiment has confirmed it.

### 3.2  Credible set and PIP — Posterior Inclusion Probability

When a GWAS hit is fine-mapped, the result is a **credible set**: a small
group of SNPs that together carry ~95% of the statistical signal at the
locus. The **PIP** of each variant is the probability that *it* is the
causal one *within that fine-mapping run*.

> **PIP = 0.8 does NOT mean "this variant is 80% likely to be causal in the
> disease". It means 80% of the posterior probability in this credible set
> sits on this variant.**

Different fine-mapping runs of the same locus (different studies, different
populations, different priors) can produce very different PIPs for the same
SNP. This is why each SNP card now shows the specific `study_id` — so you
can trace which credibleSet the PIP came from.

A locus with one credible-set SNP at PIP 0.9 is very different from a
locus with twenty SNPs each at PIP 0.05.

### 3.3  DE-z and the choice of reference cell state

The colour you see on each seed gene is its **DE-z score** — a differential
expression Z-score *relative to a reference cell state* you pick at the top.

> **DE-z = specificity vs your chosen reference. It is NOT absolute
> expression level.**

If you change the reference from "fibroblast" to "endothelial", the same
gene may flip from "highly expressed" to "barely expressed". The "right"
reference depends on your biological question:

- "What is specific to atrial cardiomyocytes vs the rest of the atrium?" →
  reference = "all other cell types in the atrium"
- "What distinguishes diseased valve cells from healthy ones?" →
  reference = "healthy valve cells"

### 3.4  AlphaGenome variant-effect score (beta)

The 🧬 floating panel shows, for select SNPs, a deep-learning prediction
of how the variant perturbs regulatory activity (chromatin accessibility,
splicing, gene expression) in each cell context.

> **This is a model prediction, not a functional assay.**

In this preview, on-demand scoring is disabled — only **15 pre-computed
representative variants** have results. If your SNP of interest is not in
the list, the panel will say so. Absence in the cache means "not yet
scored", **not** "predicted to have no effect".

### 3.5  Chat grounding and BYOK

"Grounding" means the prompt you send to the LLM is augmented with
**verified PMIDs** the backend prefetches (and, optionally, PubMed abstract
summaries the agent has previously generated and you have curated by
voting). This narrows hallucination but does not eliminate it.

**BYOK** (Bring Your Own Key) means CardioNav itself does not pay for or
store any LLM credentials. Your key is in your browser. The server uses it
transiently for one request and forgets it. See [section 6.5](#api-key-protection).

---

## 4. Interface reference

*(Detailed walkthrough of every panel and control — coming next.)*

For now, the Quickstart (section 2) covers the basics. The remaining
panels worth knowing:

- **Cytoscape graph controls** (top-right of the canvas): top-N seeds,
  peaks-per-seed, TFs-per-peak, FDR threshold sliders
- **Edge cell-type filter**: restrict edges to specific supporting cell types
- **🧬 AlphaGenome floating panel**: drag by header to reposition
- **📚 Cached refs**: backend-only optimisation (no UI exposure)

---

## 5. Chat agent

*(Coming next — the BYOK setup is in [section 2.5](#25--set-up-the-chat-agent).)*

Key points for the impatient:

- 🟢 **high** confidence = ≥ 2 net up-votes from human reviewers.
  Recurrence alone (the agent re-citing the same paper) cannot promote a
  PMID to high.
- 🟡 **medium** = 1 net vote OR ≥ 2 agent recurrences.
- ⚪ **low** = no human signal yet.
- 🔴 **flagged** = net down-votes. Not used for grounding.
- 📝 **Summarize**: explicit user gesture. Spends ~1500 of your BYOK
  tokens; saves the resulting summary as durable prior knowledge.

> ⚠️ **"Verified" badge means the PMID exists in PubMed and was linked by
> the LLM. It does NOT mean a human has read the paper and confirmed it
> supports the specific claim. Always check citations before quoting.**

---

## 6. Worked tasks

*(Five worked recipes — coming next.)*

Brief preview:

1. Find candidate TFs for the PALMD locus in AVS valve cells.
2. Pivot a CAD SNP across to AF / AVS to check overlap.
3. Compare a gene's expression across reference cell states.
4. Use the chat agent to build a literature digest for a TF–disease link.
5. Inspect AlphaGenome scores for a representative SNP.

---

## 7. FAQ & troubleshooting

*(Coming next.)*

Most common issues so far:

- **"Loading DE-z matrix" hangs for a minute or more** — likely network
  latency. From Tokyo to the Hetzner DE server the RTT is ~520 ms; from
  Cambridge UK it's ~25 ms. The DE-z files are now gzipped to ~5 MB.
- **Chat returns "⏳ Gemini key rate-limited"** — free-tier RPM/TPM is
  tight. Wait 60 s and retry, or upgrade your Google AI plan. Switching
  to a different model is also worth trying (see top-right model picker).
- **Chat returns "☁️ Gemini service overloaded (503)"** — Google capacity
  problem, not your quota. Wait 30 s and retry, your quota is intact.
- **"Summarize failed: summary not saved"** — usually the LLM returned
  malformed JSON. Just click 📝 Summarize again; the prompt has been
  hardened to reject the common Gemini failure modes.

---

## 8. Caveats & limitations

A consolidated list of the biases and traps that can lead to wrong
inferences:

- **"Verified" citation ≠ supports the claim** — see [section 5](#5-chat-agent).
- **PIP is per fine-mapping run, not per locus** — see [section 3.2](#32--credible-set-and-pip--posterior-inclusion-probability).
- **AlphaGenome covers 15 variants only** — see [section 3.4](#34--alphagenome-variant-effect-score-beta).
- **Cached refs are NOT user-curated** — the backend silently accumulates
  agent-found PMIDs to improve future grounding; this is a performance
  store, not a vetted bibliography. Only the explicit 👍 vote + 📝
  Summarize flow is "user-curated".
- **eGRN edges are predicted, not experimentally confirmed** — see
  [section 3.1](#31--egrn--enhancer-gene-regulatory-network).
- **Data is a snapshot** — OT v26.03 + GWAS sumstats current at pipeline
  build time. Variants discovered or re-fine-mapped since will not appear.

---

## 9. Glossary

| Term | One-line definition |
|---|---|
| **eGRN** | Cell-state-specific predicted regulatory network (TF → peak → gene) |
| **Seed gene** | Top GWAS-prioritised gene for a disease, shown on the outer ring |
| **ATAC peak** | Open chromatin region identified by scATAC-seq |
| **TF** | Transcription factor; predicted from sequence motif + ATAC |
| **PIP** | Posterior Inclusion Probability — credible-set membership weight |
| **L2G** | Open Targets "Locus to Gene" probabilistic SNP → gene mapping |
| **VEP** | Variant Effect Predictor (Ensembl) — most severe consequence |
| **DE-z** | Differential expression Z-score relative to a reference cell state |
| **BYOK** | Bring Your Own Key — user-supplied LLM API credentials |
| **snp2cell** | Upstream pipeline that built the per-cell-state eGRNs |
| **SCENIC+** | TF-target inference engine used inside snp2cell |
| **OT** | Open Targets — the source of credible-set + L2G annotations |
| **credibleSet** | A small group of SNPs covering ~95% of a GWAS signal at a locus |
| **study_id** | OT study identifier disambiguating which fine-mapping run produced a PIP |

---

<a name="api-key-protection"></a>
## How API keys are protected — quick reference

User-supplied Anthropic / Gemini API keys are stored **only in your
browser's `localStorage`** and sent per-request over HTTPS in an
`X-API-Key` header. The server uses them transiently (never logs,
persists, or echoes them) and scrubs any matching pattern from SDK
exception messages before returning errors. Even if the server is
compromised, no keys are at rest — they exist only in the originating
browser plus the lifetime of one in-flight HTTPS request.

---

*Questions, corrections, or feature requests → contact the PI.*
