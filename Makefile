init:
	pip install pipenv
	pip install --dev

test:
	pipenv run py.test --pep8 tests

acceptance: libs/seq_file/bin/dnacat
	$(MAKE) -C tests/from-mccortex/build

clean:
	@echo nothing to clean

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance clean
