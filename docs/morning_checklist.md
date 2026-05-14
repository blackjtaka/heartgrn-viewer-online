# Morning checklist — what to do when you wake up

Everything except the two interactive deploy steps has been pushed to GitHub. The repo at <https://github.com/blackjtaka/heartgrn-viewer-online> contains:

- Backend: `server/hub_backend.py` (Anthropic SDK, BYOK, per-IP rate limit), Dockerfile, fly.toml
- Frontend: `hub/index.html` (password gate + BYOK modal), `hub/assets/config.js`, all 195 payload JSONs
- Docs: `docs/architecture.md`, `docs/deployment.md`, this file

## Default staging credentials

- **Password (page gate)**: `heartgrn-staging-2026`
- To change it: compute `echo -n 'newpass' | shasum -a 256` and update `hub/assets/config.js` `PASSWORD_SHA256` field, then push.
- To disable the gate entirely: set `REQUIRE_PASSWORD: false` in `config.js`.

## Two steps that need you (≈10 min total)

### 1. Backend → Fly.io deploy

```bash
cd ~/Desktop/project/heartgrn-viewer-online
fly auth login          # if not already logged in
fly launch --copy-config --no-deploy   # accepts the existing fly.toml
fly deploy
```

`fly launch` will ask whether to create a Postgres / Redis — answer **No** to both. Region should default to `lhr` (London) per `fly.toml`.

Once deployed, note the URL (likely `https://heartgrn-api.fly.dev`) and verify:

```bash
curl https://heartgrn-api.fly.dev/healthz
# → {"status":"ok",...}
```

### 2. Frontend → Cloudflare Pages

Already connected via GitHub per your setup. The first push has triggered an auto-deploy. Find the URL in the Cloudflare dashboard (`https://heartgrn-viewer-online.pages.dev` or similar) and **set the build settings**:

- Framework preset: None
- Build command: (leave empty)
- Build output directory: `hub`

Then either trigger a re-deploy from the dashboard, or just push the next commit and it will re-build.

### 3. Wire frontend → backend CORS

Once both URLs are known:

```bash
# Replace with your actual URLs
fly secrets set CF_PAGES_ORIGIN=https://heartgrn-viewer-online.pages.dev
fly deploy
```

And in `hub/assets/config.js` confirm `API_BASE` matches your Fly URL (the default `heartgrn-api.fly.dev` should match if your fly app name is `heartgrn-api`). Push if changed.

## Verifying

1. Open the Cloudflare Pages URL → password prompt → enter `heartgrn-staging-2026`
2. Cytoscape graph should load (no key needed)
3. Click 🔑 Activate Agent top-right → paste a real Anthropic key → try a chat message
4. Backend `fly logs` should show the request without ever logging the key

## If something goes wrong

- **CORS errors in browser console** → the Cloudflare Pages URL is not in `ALLOWED_ORIGINS`. Set `CF_PAGES_ORIGIN` secret in Fly and redeploy. Or for a quick test, `fly secrets set ALLOW_ANY_ORIGIN=1`.
- **401 from /chat with valid key** → the SDK rejected it. Check `fly logs` for the error detail.
- **Payload 404s** → Cloudflare Pages didn't pick up `hub/` as build output. Check the build settings.
- **Frontend won't load** → look at the `pass-gate` and make sure SHA-256 is correctly set in `config.js`; the in-browser `sha256(input)` must match.

## What's left after live

- Phase 5 (Cloudflare Pages connect): you finish in dashboard
- Optional: custom domain (later)
- Optional: convert password gate → Cloudflare Access (real SSO) when ready to share with broader audience
- Optional: precompute AG1 cache for more variants and rsync into `hub/ag1_cache/` before publication

Sleep well 🌙
