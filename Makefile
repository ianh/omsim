.PHONY: clean

CFLAGS=-O2 -std=c11 -pedantic -Wall -Wno-missing-braces
LDLIBS=-lm
LLVMCC=/opt/homebrew/Cellar/llvm/12.0.1/bin/clang
EMFLAGS=--no-entry -s ALLOW_MEMORY_GROWTH=1 -s ALLOW_TABLE_GROWTH=1 --profiling-funcs -DNDEBUG
EMEXPORTS1=_malloc,_free,_verifier_create_from_bytes,_verifier_create_from_bytes_without_copying
EMEXPORTS2=_verifier_destroy,_verifier_set_cycle_limit,_verifier_error,_verifier_error_clear,_verifier_evaluate_metric
EMEXPORTS3=_verifier_set_fails_on_wrong_output,_verifier_set_fails_on_wrong_output_bonds,_verifier_wrong_output_index
EMEXPORTS4=_verifier_wrong_output_atom,_verifier_wrong_output_clear,_verifier_number_of_output_intervals
EMEXPORTS5=_verifier_output_interval,_verifier_output_intervals_repeat_after

HEADER=collision.h decode.h parse.h sim.h verifier.h
SOURCE=collision.c decode.c parse.c sim.c verifier.c

omsim: $(HEADER) $(SOURCE) Makefile main.c
	$(CC) $(CFLAGS) -o $@ $(SOURCE) main.c $(LDLIBS)

libverify.so: $(HEADER) $(SOURCE) Makefile
	$(CC) $(CFLAGS) -shared -fpic -o $@ $(SOURCE) $(LDLIBS)

libverify.wasm: $(HEADER) $(SOURCE) Makefile
	emcc $(CFLAGS) $(EMFLAGS) -s EXPORTED_FUNCTIONS=$(EMEXPORTS1),$(EMEXPORTS2),$(EMEXPORTS3),$(EMEXPORTS4),$(EMEXPORTS5) -o $@ $(SOURCE)

run-tests: $(HEADER) $(SOURCE) Makefile run-tests.c
	$(CC) $(CFLAGS) -g -D_DEFAULT_SOURCE -o $@ $(SOURCE) run-tests.c $(LDLIBS)

llvm-fuzz: $(HEADER) $(SOURCE) Makefile llvm-fuzz.c
	$(LLVMCC) $(CFLAGS) -g -fsanitize=fuzzer,address -o $@ llvm-fuzz.c $(SOURCE) $(LDLIBS)

clean:
	-rm omsim libverify.so libverify.wasm run-tests llvm-fuzz
