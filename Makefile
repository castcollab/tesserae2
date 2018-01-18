HYPOTHESIS_PROFILE = dev
PIP = pip3
BIN_DIR = "./bin"

RUN_IN_ENV = pipenv run
PYTHON = $(RUN_IN_ENV) python

UNIT_TEST_COMMAND = $(RUN_IN_ENV) pytest \
	--flake8 \
	--cov=cortexpy \
	--cov-report term-missing \
	--hypothesis-profile $(HYPOTHESIS_PROFILE)
TEST_COMMAND = BIN_DIR=$(BIN_DIR) $(UNIT_TEST_COMMAND)
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

unit:
	$(UNIT_TEST_COMMAND) cortexpy/test/test_unit

acceptance:
	$(TEST_COMMAND) cortexpy/test/test_acceptance

test:
	$(TEST_COMMAND) cortexpy

lint:
	- $(RUN_IN_ENV) pylint cortexpy \
	--ignore test \
	--disable missing-docstring,unsubscriptable-object,no-member

acceptance_: libs/seq_file/bin/dnacat
	$(MAKE) -C cortexpy/test/from-mccortex/build

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

check_git_dirty:
	git status --porcelain
	test -z "$$(git status --porcelain)"

deploy: check_git_dirty test build
	$(RUN_IN_ENV) twine upload dist/*

deploy-test: check_git_dirty test build
	$(RUN_IN_ENV) twine --repository testpypi upload dist/*

build: clean
	pipenv lock
	$(MAKE) dist

dist:
	$(PYTHON) setup.py sdist
	$(PYTHON) setup.py bdist_wheel

doc:
	$(MAKE) -C doc html

setup-benchmark:
	(test -d $(BENCHMARK_DIR) && cd $(BENCHMARK_DIR) && git pull) || git clone https://github.com/winni2k/cortex_tools_benchmark.git
	$(MAKE) -C $(BENCHMARK_DIR) setup

benchmark:
	$(MAKE) -C $(BENCHMARK_DIR)


clean:
	rm -rf dist

.PHONY: test acceptance unit clean update pipenv compile build publish doc dist
