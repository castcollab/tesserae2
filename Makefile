
test: acceptance

acceptance: libs/seq_file/bin/dnacat
	$(MAKE) -C tests/from-mccortex/build

clean:
	@echo nothing to clean

libs/seq_file/bin/dnacat:
	$(MAKE) -C $$(dirname $$(dirname $@))

.PHONY: test acceptance clean
