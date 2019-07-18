init:
	pip install pip -U
	curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
	poetry self:update
	poetry install
	poetry run pre-commit install
	poetry run pyblast install
	poetry run pyblast status


clean:
	rm -rf dist
	rm -rf pip-wheel-metadata
	rm -rf docs
	rm -rf .pytest_cache


test:
	poetry run python -m pytest


lint:
	poetry run pylint -E pydent


docs:
	echo "No documentation"


format:
	poetry run keats version up
	poetry run keats changelog up
	poetry run black shoestring tests



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
