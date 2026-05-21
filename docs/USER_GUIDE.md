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
- **65 cell states per disease** — atrial cardiomyocytes,
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

### 4.1  Top selectors and search

| Control | What it does |
|---|---|
| **Disease** | Switches the eGRN to AF / CAD / AVS. Reloads the DE-z matrix on first visit |
| **Target cell-state** | The cell state the network is *built for* — colours and edges represent regulation in this context |
| **Reference cell-state** | The denominator for DE-z. Switching this re-colours seed genes; does not rebuild the graph |
| **Score key** | Which underlying score drives the seed-gene ring (default: GWAS-z, log10 scaled) |
| **🔍 Search bar** | Cross-payload lookup: type a gene symbol or rsid → see every (disease, cs) combo where it appears |

### 4.2  Graph canvas (Cytoscape)

Concentric ring layout. From outside in:

- **Seed gene ring** — circular nodes, colour ∝ DE-z (red = up in current cs vs reference, blue = down)
- **Peak ring** — rectangular nodes labelled `chrN:start-end`. Colour ∝ GWAS-z within the peak
- **TF ring** — diamond nodes, sized by specificity score

Edges:
- **TF → peak** (inner → middle): predicted binding from SCENIC+
- **peak → gene** (middle → outer): scATAC-seq peak-to-gene linkage
- Edge thickness ∝ number of supporting cell types

Top-right sliders (above the canvas):
- **top-N seeds**: cap on outer ring count
- **peaks-per-seed**: cap on middle-ring peaks per gene
- **TFs-per-peak**: cap on inner-ring TFs per peak
- **FDR threshold**: filter low-confidence edges

> **Tip**: starting values come from the payload's `defaults`. If the graph
> looks empty or too dense, the sliders are the first thing to adjust.

### 4.3  SNP tooltip card

Click any peak / SNP node → tooltip on the right opens. Contains:

- **rsid** (clickable → dbSNP) and full **variant ID** (`chr_pos_ref_alt`)
- **PIP** with the credible-set `study_id` (Open Targets fine-mapping run)
- **gnomAD allele frequency** (joint, all populations)
- **trait** label from the GWAS study
- **VEP most-severe consequence** (e.g. `intron_variant`, `missense_variant`)
- **L2G top gene** with the OT confidence score
- **🧬 AlphaGenome variant effect** button (only enabled for the 15 beta cache variants)

### 4.4  Chat panel (right side)

- **💬 Ask Agent** toggle button at top of the side panel
- **Input** at the bottom (Cmd/Ctrl + Enter to send)
- **History** scrolls in the middle
- Citation pills appear inline; vote 👍/👎 attaches to the PMID server-side

The chat sends the current viewer state (disease, cs, visible nodes) as
context, so questions like "explain the top TF" know which TF you mean.

### 4.5  🧬 AlphaGenome floating panel

Opens when you click the 🧬 button on a SNP card. Floats above the
canvas; drag by the header to reposition. Shows:

- Pre-computed regulatory perturbation scores per track (chromatin
  accessibility, splicing, expression)
- Inline plot of the score across the surrounding ±500 kb window

**Beta**: only 15 variants are pre-computed. If you click 🧬 on a SNP
outside this set, you see the list of available variants instead.

---

## 5. Chat agent

### 5.1  BYOK key setup

Click 🔑 in the top-right toolbar. Paste either:

- **Anthropic** (`sk-ant-…`) — better instruction-following, but free
  tier doesn't exist for Anthropic; you need a paid plan
- **Google Gemini** (`AIza…`) — generous free tier:
  - `gemini-2.0-flash` (default): 15 RPM, 1M tokens/min, 1500 RPD
  - `gemini-2.5-flash` (fallback): 10 RPM, 250k TPM, 250 RPD — slightly
    higher quality, but tighter limits and more frequent 503s during
    peak hours

Where the key lives: only in your browser's `localStorage`. Sent per
request as an `X-API-Key` HTTPS header; never logged or persisted
server-side. See [How API keys are protected](#api-key-protection).

### 5.2  What the agent can do

Three broad categories of question:

1. **Lookup** — "What's the top L2G gene at the PALMD locus?"
2. **Mechanism / interpretation** — "Why might this TF matter in AVS?"
3. **View navigation** — "switch to AF, atrial cardiomyocytes, top 30 seeds"
   (the agent emits an action the viewer executes immediately)

It cannot do bulk analyses, statistical tests, or pipeline reruns.

### 5.3  Reading citations

Each cited PMID appears as a card with several pills:

| Pill | Meaning |
|---|---|
| ✓ **verified** | PMID exists on PubMed (esummary returned a real record) |
| ⚠️ **HALLUCINATED** | PMID didn't resolve. Don't trust the claim |
| **off-grounding** | Real PMID, but found via WebSearch, not in the pre-fetched grounding pool |
| 🟢 **high** | ≥ 2 net up-votes from human reviewers in this app |
| 🟡 **medium** | 1 net vote OR ≥ 2 agent recurrences in this context |
| ⚪ **low** | No human signal yet |
| 🔴 **flagged** | Net down-votes — excluded from future grounding |
| 📚 **cached** | Auto-injected from a previous chat in this context |
| 🧠 **curated** | Has a stored AI summary (see §5.5) |

> ⚠️ **"verified" only means the PMID exists. It does NOT mean a human
> has read the paper and confirmed it supports the specific claim. Always
> read citations before quoting.**

### 5.4  Voting (👍 / 👎)

Two buttons next to each citation:

- **👍** = "I read this paper and it supports the claim in context"
- **👎** = "I read this and it's wrong, off-topic, or misapplied"

Dedup: one vote per `(pmid, IP-hash)`. Re-voting overrides your previous
choice. Aggregate up/down counts drive the confidence pill above.

**Why your votes matter**: PMIDs that reach 🟢 high are auto-injected
into future chat grounding for the same `(disease, cell-state, gene)`
context. This makes the agent progressively smarter as the community
votes. The recurrence path (agent re-citing the same PMID) caps at
🟡 medium — agent-only signal cannot promote a PMID to high.

### 5.5  📝 Summarize (vote-curated prior knowledge)

On a citation card, click **📝 Summarize** to:

1. Fetch the abstract from PubMed
2. Ask your BYOK LLM to digest it into 3-4 sentences + 3-5 key findings
3. Save the digest to `/var/lib/heartgrn/curated_knowledge/<pmid>.json`

Cost: ~1500 of your BYOK tokens. Once curated, future chats in the same
context inject the digest as durable prior knowledge — the agent reads
the summary instead of re-fetching the abstract every time.

The button shows **🧠 View** once curated; click again to inline-expand
the summary card below the citation. No modal, no re-fetch.

---

## 6. Worked tasks

### 6.1  Find candidate TFs at a locus

*"Which TFs are driving the PALMD locus signal in AVS valve cells?"*

1. Top selectors: **disease = AVS**, **target cell-state = ValveInterstitial**
   (or whichever valve subset is relevant).
2. Search bar: type **`PALMD`** → it should appear as a seed-gene hit in
   the dropdown. Click it to focus the node.
3. Inspect the SNP card under PALMD: PIP, `study_id`, L2G top gene.
4. Follow edges inward from PALMD: peak nodes connected to PALMD →
   click each peak to see which TFs bind it (inner ring).
5. For each candidate TF, click its node → DE-z value tells you if the
   TF is expressed in this cell state. High specificity + binding +
   PIP-positive peak overlap = strong hypothesis.
6. Optional: in 💬 chat, ask *"What is known about &lt;TF&gt; in valve
   biology?"* and 👍 the relevant PMIDs.

### 6.2  Pivot a SNP across diseases

*"Does the rs880315 hit in AF also show up in CAD?"*

1. In the AF view, click the SNP node for rs880315 → note the
   variant ID.
2. Search bar: paste **`rs880315`** (or `1_10736809_T_C`).
3. The dropdown lists every (disease, cs) combo where this variant
   appears in a payload's SNP list. Click a CAD hit to switch.
4. Compare the surrounding network — is it the same locus driving CAD
   or has the credible set shifted?

### 6.3  Test how DE-z depends on reference

*"Is TBX5 truly cardiomyocyte-specific, or is the reference artificially
making it look that way?"*

1. Disease = AF, target = AtrialCardiomyocytesLeft.
2. Search **`TBX5`** → focus the node, note its DE-z colour.
3. Change the **Reference cell-state** dropdown to a different cell
   type (e.g. fibroblast, endothelial).
4. Watch the colour shift. If TBX5 stays strongly red across many
   references, the signal is robust. If it flips, the apparent
   specificity is reference-dependent.

### 6.4  Build a vote-curated literature digest

*"I want a persistent prior on the TBX5 → CASZ1 regulatory axis."*

1. Focus the relevant network view (the cell state where TBX5–CASZ1
   has the strongest edge support).
2. In 💬 chat: *"What evidence links TBX5 to CASZ1 regulation in
   cardiomyocytes?"*
3. Read the citations the agent returns. For each PMID you've actually
   read and find supportive: click **👍**.
4. Once a PMID reaches **🟢 high** (net up-votes ≥ 2, possibly from a
   second reviewer on a different IP), click **📝 Summarize**. Spend
   ~1500 BYOK tokens; the digest is saved server-side.
5. Subsequent chats in the same (disease, cs) context will auto-inject
   that digest as prior knowledge — the agent starts to know what you
   know.

### 6.5  Check AlphaGenome predictions for a representative SNP

*"What does AlphaGenome say about rs880315?"*

1. Focus the SNP node for rs880315 (search bar or click in the graph).
2. Click the **🧬 AlphaGenome variant effect** button on the SNP card.
3. The floating panel opens. Read the per-track scores (chromatin
   accessibility, splicing, expression) and the surrounding ±500 kb
   plot.
4. If the variant isn't in the 15-variant beta cache, the panel lists
   the variants that are — click any to inspect.

> Remember: the score is a deep-learning prediction, not a functional
> assay result. Treat it as a hypothesis to follow up experimentally.

---

## 7. FAQ & troubleshooting

### 7.1  Loading

**Q. The page shows "Loading DE-z matrix…" for a long time.**
The DE-z matrix is ~5 MB gzipped per disease. From Cambridge UK to the
Hetzner DE server it should download in ~1-2 seconds; from Tokyo it can
take 30-60 seconds due to ~520 ms round-trip latency (TCP throughput is
constrained by bandwidth-delay product on long fat pipes). Once cached,
subsequent loads are near-instant.

**Q. "Vote failed: not found".**
The `/citations/credible` route is proxied by nginx — if it's stuck on
nginx not having the `/citations/` location block, vote requests
fallthrough to a 404. Server-side fix: see the deploy notes.

### 7.2  Chat agent

**Q. "⏳ Anthropic key rate-limited".**
Anthropic Tier 1 plans cap input tokens at 30,000 / minute. CardioNav
chats average ~20-30k tokens including graph state. Wait 60 s, upgrade
plan, or switch to a Gemini key (`AIza…`).

**Q. "⏳ Gemini key rate-limited".**
Free-tier limits: `gemini-2.0-flash` 15 RPM / 1M TPM / 1500 RPD;
`gemini-2.5-flash` 10 RPM / 250k TPM / 250 RPD. Wait 60 s or upgrade
to a paid Google AI plan.

**Q. "☁️ Gemini service overloaded (503)".**
This is Google's side, not your quota — no tokens were consumed. Wait
30 s and retry. Frequent on `gemini-2.5-flash` during Asia daytime;
the default model is now `gemini-2.0-flash` for better stability.

**Q. "🔑 API key was rejected".**
Either the key prefix is wrong (`sk-ant-…` or `AIza…` required) or
the key was revoked. Re-enter via the 🔑 button.

**Q. "(empty stream)" with no message.**
Backend error mid-stream that didn't produce a structured event. Check
browser DevTools Network tab for the `/chat/stream` response. Usually
transient — retry.

### 7.3  Citations / Summarize

**Q. "Summarize failed: summary not saved (LLM may have returned
unparseable JSON; check server logs)".**
The BYOK LLM emitted output that didn't parse as JSON. Usually Gemini
wrapping the response in markdown code fences. Click 📝 Summarize
again — the prompt has been hardened to reject the common Gemini
failure modes and the JSON stripper now handles `` ```json … ``` ``
wrappers.

**Q. A "verified" PMID looks irrelevant to the claim.**
"Verified" only confirms the PMID exists. The LLM may have linked an
unrelated paper. Down-vote it (👎) to flag — the next chat in the same
context won't auto-inject it as grounding.

**Q. Why can't I see the citation history of previous users?**
Vote counts are aggregated and shown next to each PMID. Individual vote
history is anonymised (IP-hash) and not exposed in the UI.

### 7.4  Graph / SNP cards

**Q. My favourite gene isn't on the graph.**
The seed-gene ring shows the top N (default ~200) by GWAS-z. If your
gene didn't make the cut for this (disease, cs) combo, increase top-N
in the slider. If still absent, your gene was below the union-seed
threshold in the upstream pipeline.

**Q. The SNP card shows `study_id_candidates` but no single `study_id`.**
This SNP appears in multiple Open Targets credibleSets sharing the same
trait label (e.g. several AF GWAS). The payload can't pick one
unambiguously from the trait string alone — both / all candidates are
listed.

**Q. "AlphaGenome — variant not in beta cache".**
On-demand AlphaGenome scoring is disabled in this preview to keep
inference costs predictable. Only 15 representative variants have
pre-computed results; their list is shown in the modal.

### 7.5  Performance

**Q. The viewer is sluggish — laptop fan kicking in.**
The Cytoscape layout has 100-500 nodes; on integrated GPUs it can be
heavy. Try lowering the top-N sliders, or use a desktop browser
(Chrome / Firefox stable) rather than mobile / older Safari.

**Q. AlphaGenome panel won't open / shows blank.**
Hard refresh (Cmd+Shift+R / Ctrl+Shift+R) and try again. If still
blank, the variant-effect SVG fetch returned empty — the variant is
likely outside the beta cache.

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
