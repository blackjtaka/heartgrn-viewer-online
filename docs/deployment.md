# Deployment guide

Two-piece deployment: static frontend → **Cloudflare Pages**, Python backend → **Fly.io**.

## Prereqs

- GitHub repo: `blackjtaka/heartgrn-viewer-online` (this repo)
- `flyctl` installed (`brew install flyctl`)
- Cloudflare account (free)

## Backend → Fly.io

### One-off setup

```bash
fly auth signup     # or `fly auth login`
fly launch          # interactive; pick app name "heartgrn-api", region "lhr"
```

`fly launch` will detect the existing `fly.toml` and `server/Dockerfile`.

### Deploy

```bash
fly deploy
```

Verify:

```bash
curl https://heartgrn-api.fly.dev/healthz
# → {"status":"ok",...}
```

### Setting CORS allow-list

After Cloudflare Pages gives you a public URL (e.g. `https://heartgrn.pages.dev`):

```bash
fly secrets set CF_PAGES_ORIGIN=https://heartgrn.pages.dev
# Optional once you add a custom domain:
fly secrets set CUSTOM_ORIGIN=https://heartgrn.app
fly deploy
```

The backend reads these env vars at startup to populate `ALLOWED_ORIGINS`.

### Setting NCBI E-utilities key (optional, raises rate limit)

```bash
fly secrets set NCBI_API_KEY=<your-NCBI-API-key>
fly secrets set NCBI_EMAIL=your@email
```

These are server-side (used for PubMed look-ups) and have nothing to do with BYOK Anthropic.

## Frontend → Cloudflare Pages

1. https://dash.cloudflare.com → Pages → Create → Connect to Git
2. Select `blackjtaka/heartgrn-viewer-online`
3. Build settings:
   - Framework preset: None
   - Build command: (empty)
   - Build output directory: `hub`
4. Save & deploy.

Cloudflare auto-deploys on every push to `master`.

URL: `https://heartgrn.pages.dev` (or whatever Pages assigns).

### Wire frontend → backend

The frontend reads the backend URL from a constant at the top of `hub/assets/viewer.js`:

```js
const API_BASE = "https://heartgrn-api.fly.dev";
```

Update this after Fly.io deploys. Push the commit; Cloudflare re-deploys automatically.

## Custom domain (later)

1. Buy a domain (e.g. Cloudflare Registrar, ~$10/year)
2. Cloudflare dashboard → Pages → Custom domains → add `heartgrn.example.com`
3. Cloudflare auto-creates the DNS record and SSL cert
4. For the API: `fly certs add api.example.com` and add a CNAME record in Cloudflare DNS pointing to the Fly app
5. Update `API_BASE` in `viewer.js`, also `fly secrets set CUSTOM_ORIGIN=https://heartgrn.example.com`

No code change required besides the API_BASE constant.

## Rolling updates

- Code change (backend) → `fly deploy`
- Code change (frontend) → `git push` (Cloudflare auto-deploys)
- New payload data → `git push` (Cloudflare auto-deploys)
- Anthropic SDK update → bump `server/requirements.txt`, `fly deploy`

## Cost

- Fly.io free tier: 3 × 256 MB VMs (we use one). Free unless you scale up.
- Cloudflare Pages: free (unlimited bandwidth).
- LLM: $0 to operator (BYOK).
- Custom domain (optional): $10/year.
