# HeartGRN Viewer — online deployment architecture

## High-level diagram

```
┌──────────────────────────────────────────────────────────────────┐
│ User browser                                                     │
│  ├─ Cytoscape graph viewer (no key)                              │
│  └─ Agent chat panel (X-API-Key from localStorage)               │
└──────────────────────────────────────────────────────────────────┘
        │                                          │
        │ HTTPS static assets                      │ HTTPS POST /chat
        ▼                                          ▼
┌──────────────────────────┐         ┌──────────────────────────────┐
│ Cloudflare Pages         │         │ Fly.io: hub_backend          │
│ (free tier)              │         │ (Python 3.11, free VM)       │
│                          │         │                              │
│ Serves:                  │         │ Endpoints:                   │
│  hub/index.html          │         │  /chat   ── Anthropic SDK    │
│  hub/assets/*            │         │             tool-use loop    │
│  hub/payloads/*.json     │         │  /literature/{snp,tf,gene}/* │
│                          │         │             NCBI E-utilities │
│ Domain (initial):        │         │  /ag1/score/{variant}        │
│  *.pages.dev             │         │             cached results   │
│                          │         │  /healthz                    │
│ (custom domain later)    │         │                              │
└──────────────────────────┘         └──────────────────────────────┘
                                                  │
                                                  │ HTTPS (with user's key)
                                                  ▼
                                     ┌──────────────────────────────┐
                                     │ Anthropic API                │
                                     │ (claude-sonnet-4-6)          │
                                     └──────────────────────────────┘
```

## Key principles

1. **BYOK only** — user provides their own Anthropic API key for `/chat`
2. **Never persist user keys** — backend proxies the key for one request, discards it
3. **Graph is open** — Cytoscape exploration requires no key
4. **Citation grounding** — `/literature` uses NCBI E-utilities server-side to verify PMIDs the LLM mentions, preventing hallucinated citations
5. **Stateless backend** — no DB, only filesystem cache (literature, AG1)

## Cost model

| Layer | Provider | Cost |
|---|---|---|
| Static frontend | Cloudflare Pages | $0 (free tier, unlimited bandwidth) |
| Backend container | Fly.io | $0 (free tier, 3 × 256 MB VMs) |
| LLM API | Anthropic | $0 to operator (BYOK = user pays) |
| Domain | Cloudflare Registrar (optional) | $10/year |
| PubMed | NCBI | $0 (E-utilities free) |
| **Total** | | **$0–$10/year** |

## Security

- `X-API-Key` header (not cookie) — minimizes XSS surface
- CORS: allow only the Cloudflare Pages origin
- Rate limit: per-IP 20 req/min on `/chat`
- HTTPS enforced on both layers
- `Content-Security-Policy` header to prevent third-party script injection
- Anthropic key never logged, never written to disk, never echoed in responses

## Migration to custom domain (future)

1. Purchase domain at Cloudflare Registrar (~$10/year)
2. Cloudflare Pages: Custom Domains → add `heartgrn.your-domain` (auto SSL)
3. Fly.io: `fly certs add api.your-domain`
4. Update frontend `API_BASE` constant
5. Done — no other code change required
