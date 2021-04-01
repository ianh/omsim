.PHONY: clean

CFLAGS=-O2 -g -std=c11 -pedantic -Wall -Wno-missing-braces

EMFLAGS=--closure 1 -s ASSERTIONS=1 -s MODULARIZE=1 -s EXPORT_NAME=OMSIM -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall"]' -s FILESYSTEM=0 -s ALLOW_MEMORY_GROWTH=1 -s ALLOW_TABLE_GROWTH=1

wk9: wk9.c parse.c Makefile
	$(CC) $(CFLAGS) -o $@ wk9.c parse.c

web/elu-test-harness.js: wk9.c parse.c Makefile
	emcc $(CFLAGS) $(EMFLAGS) -s EXPORTED_FUNCTIONS='["_init","_load_solution","_test","_last_cycle","_last_capacity","_last_used","_last_removed","_last_query_u","_last_query_v"]' -o $@ wk9.c parse.c

web/visualize.js: visualize.c parse.c Makefile
	emcc $(CFLAGS) $(EMFLAGS) -s EXPORTED_FUNCTIONS='["_load_puzzle","_load_solution","_render","_command_buffer","_command_buffer_length","_reset","_step"]' -o $@ -Wno-gnu-folding-constant visualize.c parse.c

omsim: main.c parse.c Makefile
	$(CC) $(CFLAGS) -o $@ main.c parse.c

clean:
	rm omsim
	rm web/elu-test-harness.js
