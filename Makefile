.PHONY: clean

CFLAGS=-O2 -std=c11 -pedantic -Wall -Wno-missing-braces
LDLIBS=-lm
LLVMCC=/opt/homebrew/Cellar/llvm/12.0.1/bin/clang
EMFLAGS=--no-entry -s ALLOW_MEMORY_GROWTH=1 -s ALLOW_TABLE_GROWTH=1 -Oz --profiling-funcs -DNDEBUG
EMEXPORTS1=_malloc,_free,_verifier_create_from_bytes,_verifier_create_from_bytes_without_copying
EMEXPORTS2=_verifier_destroy,_verifier_set_cycle_limit,_verifier_error,_verifier_error_clear,_verifier_evaluate_metric
EMEXPORTS3=_verifier_set_fails_on_wrong_output,_verifier_wrong_output_index,_verifier_wrong_output_atom,_verifier_wrong_output_clear

omsim: main.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ main.c parse.c sim.c decode.c $(LDLIBS)

libverify.so: verifier.c verifier.h parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -shared -fpic -o $@ verifier.c sim.c parse.c decode.c $(LDLIBS)

libverify.wasm: verifier.c verifier.h parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	emcc $(CFLAGS) $(EMFLAGS) -s EXPORTED_FUNCTIONS=$(EMEXPORTS1),$(EMEXPORTS2),$(EMEXPORTS3) -o $@ verifier.c sim.c parse.c decode.c

run-tests: run-tests.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ run-tests.c parse.c sim.c decode.c $(LDLIBS)

llvm-fuzz: llvm-fuzz.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(LLVMCC) $(CFLAGS) -fsanitize=fuzzer,address -o $@ llvm-fuzz.c parse.c sim.c decode.c $(LDLIBS)

clean:
	-rm omsim libverify.so libverify.wasm run-tests llvm-fuzz
