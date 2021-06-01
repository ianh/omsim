.PHONY: clean

CFLAGS=-O2 -std=c11 -pedantic -Wall -Wno-missing-braces -g
LDLIBS=-lm

omsim: main.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ main.c parse.c sim.c decode.c $(LDLIBS)

libverify.so: verifier.c verifier.h parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -shared -fpic -o $@ verifier.c sim.c parse.c decode.c $(LDLIBS)

run-tests: run-tests.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ run-tests.c parse.c sim.c decode.c $(LDLIBS)

clean:
	rm omsim
