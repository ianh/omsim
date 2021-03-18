.PHONY: clean

CFLAGS=-O2 -g -std=c11 -pedantic -Wall -Wno-missing-braces

omsim: main.c parse.c Makefile
	$(CC) $(CFLAGS) -o $@ main.c parse.c

clean:
	rm omsim
