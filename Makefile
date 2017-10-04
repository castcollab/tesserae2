HYPOTHESIS_PROFILE = dev

init: update
	pip install pipenv
	pipenv install --dev
	$(MAKE) -C libs/mccortex all MAXK=63
	$(MAKE) -C libs/mccortex all MAXK=31

update:
	git submodule update --init --recursive

test:
	pipenv run pytest pycortex \
	--flake8 \
	--cov=pycortex \
	--cov-report term-missing \
	--hypothesis-profile $(HYPOTHESIS_PROFILE)

acceptance: libs/seq_file/bin/dnacat
	$(MAKE) -C pycortex/test/from-mccortex/build

clean:
	@echo nothing to clean

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance clean
