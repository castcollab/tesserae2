HYPOTHESIS_PROFILE = dev
PIP = pip3
BIN_DIR := "./bin"

RUN_IN_ENV = pipenv run
PYTHON = $(RUN_IN_ENV) python

FAST_TEST_COMMAND = $(RUN_IN_ENV) pytest \
	--hypothesis-profile dev

PYTEST_COMMAND = pytest \
           --flake8 \
           --cov cortexpy \
           --cov-report term-missing \
           --cov-report html \
           --cov-config setup.cfg \
           --hypothesis-profile $(HYPOTHESIS_PROFILE)
BASE_TEST_COMMAND = PATH=$(BIN_DIR):$(PATH) tox -- $(PYTEST_COMMAND)
TEST_COMMAND = $(BASE_TEST_COMMAND) tests/test_unit tests/test_acceptance
BENCHMARK_DIR := cortex_tools_benchmark

init:
	$(MAKE) pipenv
	$(MAKE) compile

dev: init
	pipenv install -e .
	echo development environment set up at $$(pipenv --venv)

pipenv: update
	$(PIP) install pipenv
	pipenv install --dev

compile: update
	$(MAKE) -C libs/mccortex all MAXK=63
	$(MAKE) -C libs/mccortex all MAXK=31

update:
	git submodule update --init --recursive

doctest:
	tox -- python -m doctest $$(find src/cortexpy/graph -type f -name '*.py' | grep -v '__*' )

fast:
	$(FAST_TEST_COMMAND) tests/test_unit

unit:
	$(BASE_TEST_COMMAND) tests/test_unit

acceptance:
	$(BASE_TEST_COMMAND) tests/test_acceptance

fixtures:
	$(MAKE) -C $(BENCHMARK_DIR) test-fixtures

pycompile:
	$(PYTHON) setup.py build_ext --inplace

test:
	$(MAKE) doctest
	$(TEST_COMMAND)

ci:
	$(TEST_COMMAND) -vv

lint:
	tox -- pytest --pylint src

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

check_git_dirty:
	git status --porcelain
	test -z "$$(git status --porcelain)"

deploy: check_git_dirty
	tox -r
	$(MAKE) test
	$(MAKE) build
	$(MAKE) upload

upload:
	$(RUN_IN_ENV) twine upload dist/*

build: clean
	$(MAKE) check_git_dirty
	$(MAKE) dist

dist: pycompile
	$(PYTHON) setup.py sdist

doc:
	$(MAKE) -C doc html

setup-benchmark: dist
	$(eval CORTEXPY_WHEEL := $(shell find dist/cortexpy-*.whl))
	CORTEXPY_WHEEL=../$(CORTEXPY_WHEEL) $(MAKE) -C $(BENCHMARK_DIR) setup

benchmark:
	$(MAKE) -C $(BENCHMARK_DIR)

clean:
	rm -rf dist

.PHONY: test acceptance unit clean update pipenv compile build publish doc dist pycompile
