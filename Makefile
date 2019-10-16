PIP=pip3

.PHONY: docs  # necessary so it doesn't look for 'docs/makefile html'

init:
	pip install pip -U
	curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
	poetry self:update
	poetry install
	poetry run pre-commit install
	poetry run pyblast install
	poetry run pyblast status


check:
	poetry run pre-commit run


clean:
	rm -rf dist
	rm -rf pip-wheel-metadata
	rm -rf docs
	rm -rf .pytest_cache


test:
	poetry run python -m pytest


lint:
	poetry run pylint -E pydent


pullversion:
	poetry run keats version up


docs: | pullversion
	@echo "Updating docs"

	# copy README.md to README.rst format for Sphinx documentation
	# we can comment this out if we do not want to include the README.md in the sphinx documentation

	poetry run keats changelog up
	cp .keats/changelog.md docsrc/developer/changelog.md

	rm -rf docs
	cd docsrc && poetry run make html
	find docs -type f -exec chmod 444 {} \;
	@echo "\033[95m\n\nBuild successful! View the docs homepage at docs/html/index.html.\n\033[0m"

	touch docs/.nojekyll
	open ./docs/index.html


format:
	poetry run keats version up
	poetry run keats changelog up
	poetry run black dasi tests



lock:
	poetry run upver
	poetry update


build:
	poetry run upver
	poetry build


release:
	sh scripts/release.sh


klocs:
	find . -name '*.py' | xargs wc -l


colab:
	cd scratch; jupyter notebook \
   --NotebookApp.allow_origin='https://colab.research.google.com' \
   --NotebookApp.port_retries=0 --port=8889


updatedeps:
	poetry cache:clear --all pypi
	poetry update


benchmark:
	poetry run pytest tests/test_benchmark --benchmark-autosave --benchmark-max-time=0.1
	poetry run pytest-benchmark compare
