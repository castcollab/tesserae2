HYPOTHESIS_PROFILE = dev

init: update
	pip install pipenv
	pipenv install --dev
	$(MAKE) -C libs/mccortex all MAXK=63
	$(MAKE) -C libs/mccortex all MAXK=31

update:
	git submodule update --init --recursive

unit:
	pipenv run pytest pycortex/test/test_unit \
	--hypothesis-profile $(HYPOTHESIS_PROFILE)

test:
	pipenv run pytest pycortex \
	--flake8 \
	--cov=pycortex \
	--cov-report term-missing \
	--hypothesis-profile $(HYPOTHESIS_PROFILE)

lint:
	- pipenv run pylint pycortex \
	--ignore test \
	--disable missing-docstring,unsubscriptable-object,no-member

acceptance: libs/seq_file/bin/dnacat
	$(MAKE) -C pycortex/test/from-mccortex/build

clean:
	@echo nothing to clean

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance clean
