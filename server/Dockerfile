# Backend image for hub_backend (BYOK Anthropic + PubMed proxy)
FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# Install minimal system deps (none beyond stdlib for now; keep slim)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Python deps
COPY server/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# App code + static hub content
COPY server/hub_backend.py /app/server/hub_backend.py
COPY hub /app/hub

# Non-root user
RUN useradd --uid 1001 --create-home --shell /bin/false appuser \
 && chown -R appuser:appuser /app
USER appuser

EXPOSE 8766

# Bind to 0.0.0.0 so Fly.io can reach the listener
CMD ["python", "-u", "server/hub_backend.py", \
     "--bind", "0.0.0.0", \
     "--port", "8766", \
     "--hub-dir", "/app/hub"]
