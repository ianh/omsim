.PHONY: clean

CFLAGS=-O2 -std=c11 -pedantic -Wall -Wno-missing-braces -g
LDLIBS=-lm
LLVMCC=/opt/homebrew/Cellar/llvm/12.0.0/bin/clang
EMFLAGS=--no-entry -s ALLOW_MEMORY_GROWTH=1 -s ALLOW_TABLE_GROWTH=1

omsim: main.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ main.c parse.c sim.c decode.c $(LDLIBS)

libverify.so: verifier.c verifier.h parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -shared -fpic -o $@ verifier.c sim.c parse.c decode.c $(LDLIBS)

libverify.wasm: verifier.c verifier.h parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	emcc $(CFLAGS) $(EMFLAGS) -o $@ verifier.c sim.c parse.c decode.c

run-tests: run-tests.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ run-tests.c parse.c sim.c decode.c $(LDLIBS)

llvm-fuzz: llvm-fuzz.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(LLVMCC) $(CFLAGS) -fsanitize=fuzzer,address -o $@ llvm-fuzz.c parse.c sim.c decode.c $(LDLIBS)

clean:
	rm omsim
	rm libverify.so
	rm libverify.wasm
	rm run-tests
	rm llvm-fuzz
