init:
	pip install pipenv
	pipenv install --dev

test:
	pipenv run py.test --pep8 pycortex

acceptance: libs/seq_file/bin/dnacat
	$(MAKE) -C pycortex/test/from-mccortex/build

clean:
	@echo nothing to clean

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance clean
