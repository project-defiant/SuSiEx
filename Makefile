.PHONY: help dev lint test native-test clean

help: ## Show available targets
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-12s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

dev: ## Install development dependencies
	@uv sync --extra test

lint: ## Run linting, formatting, and type checks
	@uv run ruff check src/susiex_cli tests/python
	@uv run ruff format --check src/susiex_cli tests/python
	@uv run ty check src/susiex_cli

test: ## Run the Python test suite
	@uv run pytest -q

native-test: ## Run the C++ regression suite
	@bash run_tests.sh

clean: ## Remove local Python build and test artifacts
	@rm -rf .venv .pytest_cache .ruff_cache build dist *.egg-info
