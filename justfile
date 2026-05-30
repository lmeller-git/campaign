set dotenv-load := true
DOCKER_USER := "lmeller"

HF_TOKEN_VAL := env_var_or_default("HF_TOKEN", "")
HF_HOME := env_var_or_default("HF_HOME", "")
HF_HUB_CACHE := env_var_or_default("HF_HUB_CACHE", "")
DOCKER_PASSWORD_VAL := env_var_or_default("DOCKER_PASSWORD", "")

default:
    @just --list

auth:
    uvx --from lyceum-cli lyceum auth login

setup:
    git submodule update --init --remote --recursive --force
    uv sync --frozen

run *args:
    uv run python python/main.py {{args}}

build:
    @export HF_TOKEN=$(grep HF_TOKEN .env | cut -d '"' -f 2) && \
    docker build --secret id=HF_TOKEN,env=HF_TOKEN -t lmeller/campaign-img:latest .

run-container-cpu:
    docker run --rm -it --env-file .env -v ./data:/app/data lmeller/campaign-img:latest

run-container-gpu:
    docker run --gpus all --rm -it --env-file .env -v ./data:/app/data lmeller/campaign-img:latest

lint:
    uv run ruff format .
    uv run ruff check .

test:
    uv run pytest

requirements: setup
    uv export --no-emit-workspace --no-dev --no-annotate --no-header --no-hashes --output-file requirements.txt

fetch_result execution_id:
    curl https://api.lyceum.technology/api/v2/external/execution/{{execution_id}} -H "Authorization: Bearer $LYCEUM_API_KEY"

fetch_log execution_id:
    curl https://api.lyceum.technology/api/v2/external/logs/execution/{{execution_id}} -H "Authorization: Bearer $LYCEUM_API_KEY"

run-cloud machine="cpu" *args="": requirements
    uvx --from lyceum-cli lyceum python run -m {{machine}} -r requirements.txt python/main.py {{args}}


run-cloud-image machine="gpu.h100":
    uvx --from lyceum-cli lyceum docker run {{DOCKER_USER}}/campaign-img:latest \
      -m {{machine}} \
      -e HF_TOKEN="{{HF_TOKEN_VAL}}" \
      -e HF_HOME="{{HF_HOME}}" \
      -e HF_HUB_CACHE="{{HF_HUB_CACHE}}" \
      --registry-type basic \
      --registry-creds '{"username":"{{DOCKER_USER}}","password":"{{DOCKER_PASSWORD_VAL}}"}'

cloud-logs exec_id:
    uvx --from lyceum-cli lyceum docker logs {{exec_id}}

