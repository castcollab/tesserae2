HYPOTHESIS_PROFILE = dev
PIP = pip3
BIN_DIR = "./bin"

RUN_IN_ENV = pipenv run
PYTHON = $(RUN_IN_ENV) python

FAST_TEST_COMMAND = pytest \
	--hypothesis-profile $(HYPOTHESIS_PROFILE)

BASE_TEST_COMMAND = $(RUN_IN_ENV) tox -- pytest \
           --flake8 \
           --cov \
           --cov-report term-missing \
           --cov-report html \
           --hypothesis-profile $(HYPOTHESIS_PROFILE)
TEST_COMMAND = $(BASE_TEST_COMMAND) tests/test_unit tests/test_acceptance
BENCHMARK_DIR := cortex_tools_benchmark

init: update pipenv compile

pipenv: update
	$(PIP) install pipenv
	pipenv install --dev

compile: update
	$(MAKE) -C libs/mccortex all MAXK=63
	$(MAKE) -C libs/mccortex all MAXK=31

update:
	git submodule update --init --recursive

fast:
	$(FAST_TEST_COMMAND) tests/test_unit

unit:
	$(BASE_TEST_COMMAND) tests/test_unit

acceptance:
	BIN_DIR=$(BIN_DIR) $(BASE_TEST_COMMAND) tests/test_acceptance

fixtures:
	$(MAKE) -C $(BENCHMARK_DIR) test-fixtures

pycompile:
	$(PYTHON) setup.py build_ext --inplace

test: pycompile
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
