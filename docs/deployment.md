# Deployment guide

Single-host deployment: **nginx** serves the static frontend (`hub/`) and reverse-proxies API paths to the Python backend (`server/hub_backend.py`), which runs as a **systemd** service on `127.0.0.1:8766`. See [architecture.md](architecture.md) for the overall picture.

## Prereqs

- A Linux host with nginx, certbot, and Python ≥ 3.11
- A DNS name pointing at the host (examples below use `viewer.example.org`)

## 1. Code and Python environment

```bash
sudo mkdir -p /opt/heartgrn && sudo chown "$USER" /opt/heartgrn
git clone https://github.com/blackjtaka/heartgrn-viewer-online.git /opt/heartgrn
cd /opt/heartgrn
python3 -m venv .venv
.venv/bin/pip install -r server/requirements.txt
```

nginx serves static files from `/opt/heartgrn/hub` (the `root` in `deploy/nginx.conf.sample`).

## 2. Backend → systemd

The backend binds to `127.0.0.1` by default and serves `hub/` next to `server/` (override with `--hub-dir`). The unit file is not in this repo; a minimal example:

```ini
# /etc/systemd/system/heartgrn-backend.service
[Unit]
Description=CardioNav backend (hub_backend.py)
After=network-online.target

[Service]
WorkingDirectory=/opt/heartgrn
ExecStart=/opt/heartgrn/.venv/bin/python server/hub_backend.py --port 8766
Environment=ALLOWED_ORIGIN=https://viewer.example.org
Environment=CUSTOM_ORIGIN=https://www.heartcellatlas.org
EnvironmentFile=-/etc/heartgrn/backend.env
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now heartgrn-backend
curl -s http://127.0.0.1:8766/healthz
# → {"status": "ok", ...}
```

### Environment variables

| Variable | Purpose |
|---|---|
| `ALLOWED_ORIGIN`, `CUSTOM_ORIGIN` | Origins allowed by the CORS / Origin check (localhost is always allowed) |
| `CHAT_RATE_PER_MIN` | Backend per-IP limit on `/chat` (default 20) |
| `NCBI_API_KEY`, `NCBI_EMAIL` | Optional; raise the NCBI E-utilities rate limit from 3/s to 10/s and identify the caller. Server-side only, unrelated to BYOK |
| `GEMINI_MODEL` | Default Gemini model (default `gemini-2.5-flash`) |
| `CREDIBLE_CITATIONS_FILE` | Credible-citation log (default `/var/lib/heartgrn/credible_citations.jsonl`) |
| `CREDIBLE_IP_SALT` | Salt for hashing voter IPs — set a private value in production |
| `AG1_PYTHON`, `AG1_SCRIPT` | Interpreter and `ag1.py` path for on-demand AlphaGenome scoring; only needed for variants not in the precomputed cache |

Put private values (`NCBI_API_KEY`, `CREDIBLE_IP_SALT`) in `/etc/heartgrn/backend.env`, not in the repo.

## 3. nginx

Follow the header of [`deploy/nginx.conf.sample`](../deploy/nginx.conf.sample):

1. Replace `viewer.example.org` with your domain
2. Copy to `/etc/nginx/sites-available/heartgrn` and symlink into `sites-enabled/`
3. Create Basic Auth users: `sudo htpasswd -c /etc/nginx/heartgrn.htpasswd <user>`
4. Add the chat rate-limit zone to `/etc/nginx/nginx.conf` (http block):
   `limit_req_zone $binary_remote_addr zone=chat:10m rate=30r/m;`
5. `sudo certbot --nginx -d <your-domain>`
6. `sudo nginx -t && sudo systemctl reload nginx`

`/healthz` is exempt from Basic Auth for monitoring; everything else requires it.

## 4. Serving under a path prefix

The public URL is `https://www.heartcellatlas.org/cardionav/`. The Sanger load balancer strips `/cardionav/` before forwarding, so nginx sees root paths. In the browser, `hub/assets/config.js` detects the prefix and sets `API_BASE` to `"/cardionav"`; no rebuild is needed. For local development against a separately running backend, set `API_BASE: "http://127.0.0.1:8766"` in `config.js`.

## Rolling updates

```bash
cd /opt/heartgrn && git pull
# backend code or requirements changed:
.venv/bin/pip install -r server/requirements.txt
sudo systemctl restart heartgrn-backend
# frontend / payload changes take effect on the next page load (static files)
```
