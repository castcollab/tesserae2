HYPOTHESIS_PROFILE = dev
PIP = pip

TEST_COMMAND = pipenv run pytest \
	--flake8 \
	--cov=pycortex \
	--cov-report term-missing \
	--hypothesis-profile $(HYPOTHESIS_PROFILE)


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
	$(TEST_COMMAND) pycortex/test/test_unit

acceptance:
	$(TEST_COMMAND) pycortex/test/test_acceptance

test:
	$(TEST_COMMAND) pycortex


lint:
	- pipenv run pylint pycortex \
	--ignore test \
	--disable missing-docstring,unsubscriptable-object,no-member

acceptance_: libs/seq_file/bin/dnacat
	$(MAKE) -C pycortex/test/from-mccortex/build

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance unit clean update pipenv compile
