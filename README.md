# HeartGRN Viewer — online

Public web deployment of the [HeartGRN viewer](https://github.com/blackjtaka/heartgrn-viewer), an interactive Cytoscape-based explorer for cardiac GWAS-eGRN propagation results from the snp2cell project.

## Live (work in progress)

- Frontend: TBD (Cloudflare Pages)
- Backend: TBD (Fly.io)

## What's free, what needs an API key

| Feature | Requires API key |
|---|---|
| Cytoscape graph explorer (nodes / edges / paths / panels) | **No** — fully open |
| Literature panel (PubMed E-utilities) | No (server cache) |
| Variant ATAC/RNA tracks (AG1 cache) | No (precomputed) |
| **Agent chat** (interactive LLM Q&A about the graph) | **Yes** — bring your own Anthropic API key |

## Bring your own key (BYOK)

The agent chat uses your own Anthropic API key, supplied at activation:

1. Sign up at [console.anthropic.com](https://console.anthropic.com) (new users get $5 free credit ≈ 1,000 graph queries)
2. Create an API key
3. Click "Activate Agent" on the viewer and paste your key
4. The key is stored in your browser's localStorage only — never sent to or stored on our servers

The backend only proxies the key to Anthropic on each `/chat` request and discards it after the response.

## Project structure

```
heartgrn-viewer-online/
├── hub/              ← static frontend → Cloudflare Pages
│   ├── index.html
│   ├── assets/
│   │   ├── viewer.js
│   │   └── styles.css
│   └── payloads/     ← per (disease, cell_state) JSON
├── server/           ← Python backend → Fly.io
│   ├── hub_backend.py
│   ├── requirements.txt
│   └── Dockerfile
├── docs/
│   └── architecture.md
└── scripts/
    └── sync_from_main.sh   ← rsync payloads from upstream snp2cell repo
```

## Related work

This viewer accompanies the cardiac GWAS propagation analysis in our snp2cell project. See the [upstream development repo](https://github.com/blackjtaka/heartgrn-viewer) for the analysis pipeline and figure-generation scripts.

## License

MIT — see [LICENSE](LICENSE).

## Contact

Takahiro Jimba (tj372@cam.ac.uk) — University of Cambridge.
