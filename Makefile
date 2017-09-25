HYPOTHESIS_PROFILE = dev

PIP=pip

init: update
	$(PIP) install pipenv
	pipenv install --dev
	$(MAKE) -C libs/mccortex all MAXK=63
	$(MAKE) -C libs/mccortex all MAXK=31

update:
	git submodule update --init --recursive
	cd libs/mccortex/libs/htslib && git checkout 49fdfbd # check out 1.5.0 release

test:
	pipenv run pytest pycortex --pep8 --hypothesis-profile $(HYPOTHESIS_PROFILE)

acceptance: libs/seq_file/bin/dnacat
	$(MAKE) -C pycortex/test/from-mccortex/build

clean:
	@echo nothing to clean

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance clean
