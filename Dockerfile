FROM pytorch/pytorch:2.12.0-cuda12.6-cudnn9-devel AS base

RUN apt-get update && apt-get install -y \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    curl -LsSf https://just.systems/install.sh | bash -s -- --to /usr/local/bin

WORKDIR /app

ENV PATH="/app/.venv/bin:$PATH"
ENV HF_HOME="/app/data/.hf_cache"
ENV HF_HUB_CACHE="/app/data/.hf_cache/hub"

RUN mkdir -p /app/data/.hf_cache && chmod -R 777 /app/data

COPY pyproject.toml uv.lock README.md ./

RUN --mount=type=secret,id=HF_TOKEN \
    export HF_TOKEN=$(cat /run/secrets/HF_TOKEN) && \
    uv sync --frozen --no-dev

ENTRYPOINT ["run_pipeline.sh"]
