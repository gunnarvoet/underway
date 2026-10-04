.PHONY: check format format-check docs ghdocs servedocs test help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := uv run python -c "$$BROWSER_PYSCRIPT"

help:
	@uv run python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

check: ## check style
	uv run ruff check src tests

format: ## format code using ruff
	uv run ruff format src tests

format-check: ## check code style using ruff format --diff
	uv run ruff format --diff src tests

docs: ## generate documentation using pdoc and open it
	rm -rf docs
	uv run pdoc -d numpy -o docs -t .pdoc-theme-gv --math src/underway/
	$(BROWSER) docs/index.html

ghdocs: ## generate documentation using pdoc (used by the GitHub workflow)
	rm -rf docs
	uv run pdoc -d numpy -o docs -t .pdoc-theme-gv --math src/underway/

servedocs: ## serve the docs & watch for changes
	uv run pdoc -d numpy -t .pdoc-theme-gv --math src/underway

test: ## run tests
	uv run pytest
