HYPOTHESIS_PROFILE = dev
PIP = pip3
BIN_DIR = "./bin"

UNIT_TEST_COMMAND = pipenv run pytest \
	--flake8 \
	--cov=cortexpy \
	--cov-report term-missing \
	--hypothesis-profile $(HYPOTHESIS_PROFILE)
TEST_COMMAND = BIN_DIR=$(BIN_DIR) $(UNIT_TEST_COMMAND)

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
	- pipenv run pylint cortexpy \
	--ignore test \
	--disable missing-docstring,unsubscriptable-object,no-member

acceptance_: libs/seq_file/bin/dnacat
	$(MAKE) -C cortexpy/test/from-mccortex/build

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance unit clean update pipenv compile
