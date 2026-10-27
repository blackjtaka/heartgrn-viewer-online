# CardioNav (HeartGRN Viewer) — online

Web deployment of CardioNav, an interactive Cytoscape-based explorer for cardiac GWAS → eGRN propagation results (AF / CAD / AVS × heart cell states) from the snp2cell project. See [docs/USER_GUIDE.md](docs/USER_GUIDE.md) for how to use it.

## Live

- Preview: <https://www.heartcellatlas.org/cardionav/> (access by invitation; HTTP Basic Auth)

## What's free, what needs an API key

| Feature | Requires API key |
|---|---|
| Cytoscape graph explorer (nodes / edges / paths / panels) | **No** |
| Literature panel (PubMed E-utilities) | No (server cache) |
| Variant ATAC/RNA tracks (AG1 cache) | No (precomputed) |
| **Agent chat** (interactive LLM Q&A about the graph) | **Yes** — bring your own Anthropic API key |

## Bring your own key (BYOK)

The agent chat uses your own Anthropic API key, supplied at activation:

1. Sign up at [console.anthropic.com](https://console.anthropic.com) and create an API key
2. Click "Activate Agent" on the viewer and paste your key
3. The key is kept in your browser's localStorage

On each `/chat` request the key is sent to the backend, which forwards it to the provider and discards it after the response; key-shaped strings are redacted from server error messages. A Google Gemini key (`AIza…`) is also accepted.

## Project structure

```
heartgrn-viewer-online/
├── hub/                  ← static frontend (served by nginx)
│   ├── index.html
│   ├── assets/           ← viewer.js, config.js, specificity.js, styles.css, logo
│   ├── payloads/         ← per (disease, cell_state) JSON
│   ├── manifest.json
│   ├── lineage.json      ← cell-state hierarchy + cell counts
│   └── snp_annotations.json
├── server/               ← Python backend (hub_backend.py, requirements.txt, Dockerfile)
├── deploy/
│   └── nginx.conf.sample ← reverse-proxy template (replace viewer.example.org)
├── docs/                 ← USER_GUIDE.md, architecture.md, deployment.md
└── scripts/
    └── add_study_id_to_payloads.py
```

### Backend environment variables

| Variable | Purpose |
|---|---|
| `ALLOWED_ORIGIN` / `CUSTOM_ORIGIN` | CORS origins of the deployed frontend |
| `AG1_PYTHON`, `AG1_SCRIPT` | Python interpreter and `ag1.py` path for on-demand AlphaGenome scoring (only needed for variants not in the precomputed cache) |

## Related work

This viewer accompanies the cardiac GWAS propagation analysis in our snp2cell project. The analysis pipeline and figure-generation scripts live in a separate repository that is not yet public.

## License

MIT — see [LICENSE](LICENSE).

## Contact

Takahiro Jimba (tj372@cam.ac.uk) — University of Cambridge.
