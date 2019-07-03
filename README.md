# PoetryBoilerPlate

This repo provides boilerplate code and structure for new python projects. Just fork to begin your new python package.

**Tool chain**

* `poetry` - package manager
* `black` - python formatter
* `pre-commit` - git commit hooks for autoformatting
* `pytest` - testing framework
* `make` - common commands

## Changing your package name

Change the name `pkg` to your package name at the following locations

* pyproject.toml[tools.poetry][name]
* pyproject.toml[tools.poetry][scripts]
* folder `pkg`
* `Makefile` at `format` function

## Install tools

This code uses **poetry** https://poetry.eustace.io/ to manage installation. Version is managed in the `pypoetry.toml`
file. Be sure to edit this file and be sure to change the `name` key to the name of your package.

To install your package, run:

```
make
```

## Version control

The true version of your package is maintained in the `pypoetry.toml` file. From within your package,
version is maintained in a json file in `pkg/_version/version.json`. To update this file, run:

```
poetry run version
```

## Formatting

Formatting is maintained by `black` https://black.readthedocs.io/en/stable/. A git hook is

## Interactive release to PyPI

You may run an interactive release script using:

```
make release
```

If you'd like to release to an alternative repository, look at the **poetry** documentation
to edit your configuration file:

```
poetry config repository.<REPO> <URL>
```

