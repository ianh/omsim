.PHONY: clean

CFLAGS=-O2 -std=c11 -pedantic -Wall -Wno-missing-braces -g

EMFLAGS=--closure 1 -s MODULARIZE=1 -s EXPORT_NAME=OMSIM -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall"]' -s FILESYSTEM=0 -s ALLOW_MEMORY_GROWTH=1 -s ALLOW_TABLE_GROWTH=1

wk9: wk9.c parse.c Makefile
	$(CC) $(CFLAGS) -o $@ wk9.c parse.c

libverify.so: libverify.c parse.c sim.c decode.c parse.h sim.h decode.h area.c area.h Makefile
	$(CC) $(CFLAGS) -shared -fpic -o $@ libverify.c sim.c parse.c decode.c area.c

web/elu-test-harness.js: wk9.c parse.c Makefile
	emcc $(CFLAGS) $(EMFLAGS) -s EXPORTED_FUNCTIONS='["_init","_load_solution","_test","_last_cycle","_last_capacity","_last_used","_last_removed","_last_query_u","_last_query_v"]' -o $@ wk9.c parse.c

web/draw.js: draw.c parse.c sim.c Makefile
	emcc $(CFLAGS) $(EMFLAGS) -s EXPORTED_FUNCTIONS='["_init","_load_solution","_test"]' -o $@ draw.c parse.c sim.c

web/visualize.js: visualize.c parse.c sim.c parse.h sim.h Makefile
	emcc $(CFLAGS) $(EMFLAGS) -s EXPORTED_FUNCTIONS='["_load_puzzle","_load_solution","_render","_command_buffer","_command_buffer_length","_reset","_step"]' -o $@ -Wno-gnu-folding-constant visualize.c parse.c sim.c

omsim: main.c parse.c sim.c decode.c area.c parse.h sim.h decode.h area.h Makefile
	$(CC) $(CFLAGS) -o $@ main.c parse.c sim.c decode.c area.c

run-tests: run-tests.c parse.c sim.c decode.c area.c parse.h sim.h decode.h area.h Makefile
	$(CC) $(CFLAGS) -o $@ run-tests.c parse.c sim.c decode.c area.c

clean:
	rm omsim
	rm web/elu-test-harness.js
