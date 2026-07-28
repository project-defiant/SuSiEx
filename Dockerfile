FROM python:3.11-slim

COPY --from=ghcr.io/astral-sh/uv:0.10.7 /uv /uvx /bin/

RUN apt-get update \
    && apt-get install --no-install-recommends -y build-essential procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY pyproject.toml uv.lock setup.py README.md LICENSE /app/
COPY include /app/include
COPY src /app/src

RUN uv sync --frozen --no-dev

ENV PATH="/app/.venv/bin:${PATH}" \
    PYTHONDONTWRITEBYTECODE=1

ENTRYPOINT ["susiex"]
