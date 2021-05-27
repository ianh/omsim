.PHONY: clean

CFLAGS=-O2 -std=c11 -pedantic -Wall -Wno-missing-braces -g

omsim: main.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ main.c parse.c sim.c decode.c

libverify.so: verifier.c verifier.h parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -shared -fpic -o $@ verifier.c sim.c parse.c decode.c

run-tests: run-tests.c parse.c sim.c decode.c parse.h sim.h decode.h Makefile
	$(CC) $(CFLAGS) -o $@ run-tests.c parse.c sim.c decode.c

clean:
	rm omsim
