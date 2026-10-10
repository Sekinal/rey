FROM ghcr.io/astral-sh/uv:0.11.15 AS uv
FROM debian:bookworm-slim
COPY --from=uv /uv /usr/local/bin/uv
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates texlive-latex-base texlive-latex-recommended texlive-fonts-recommended lmodern poppler-utils && rm -rf /var/lib/apt/lists/*
WORKDIR /proyecto
COPY . .
RUN uv sync --locked
CMD ["sh", "-c", "uv run --locked pytest && uv run --locked feal"]
