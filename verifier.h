#ifndef OM_VERIFIER_H
#define OM_VERIFIER_H

// this "verifier" API is designed to be called from an FFI.

// takes the path of a .puzzle and a .solution file.  returns a pointer to a
// newly-allocated verifier object.  if the files couldn't be parsed, the
// verifier_error() function will return a non-null string value describing the
// error.
void *verifier_create(const char *puzzle_filename, const char *solution_filename);

// like verifier_create(), but takes two byte arrays instead of reading from the
// filesystem.
void *verifier_create_from_bytes(const char *puzzle_bytes, int puzzle_length,
 const char *solution_bytes, int solution_length);

// like verifier_create_from_bytes(), but doesn't make a copy of the bytes.  the
// caller is responsible for keeping these allocations around until the verifier
// is destroyed.
void *verifier_create_from_bytes_without_copying(const char *puzzle_bytes, int puzzle_length,
 const char *solution_bytes, int solution_length);

// destroys a verifier object created by verifier_create(), freeing its memory.
void verifier_destroy(void *verifier);

// set how many cycles to wait for a solution to complete.
void verifier_set_cycle_limit(void *verifier, int cycle_limit);

// returns a string describing any error that occurred during solution parsing,
// decoding, or verification.  if there haven't been any errors, returns a null
// pointer.  returned strings are valid for the lifetime of the program and
// don't need to be freed.  calling this function doesn't clear the error;
// to clear the error, create a new verifier or call verifier_error_clear().
const char *verifier_error(void *verifier);

// after calling this function, verifier_error() will return null again.
void verifier_error_clear(void *verifier);

// verifier_evaluate_metric() evaluates a metric for the given solution.  if
// there is an error, verifier_evaluate_metric() returns -1 and sets an error
// string, which verifier_error() will return.
//
// the "parsed cycles", "parsed cost", "parsed area", and "parsed instructions"
//  metrics return the value recorded in the solution file.  if no value is
//  recorded, returns zero.
// the "product <product number> <metric>" metric measures <metric> as if the
//  solution completes after the <product number>-th product instead of after
//  the product count specified in the puzzle file.  for example,
//  "product 1 cycles" will measure the latency for the first product.
// the "parts of type <part type>" metrics return the number of parts in the
//  solution file which have that type.  here's a list of the possible types:
//   "parts of type arm1"
//   "parts of type arm2"
//   "parts of type arm3"
//   "parts of type arm6"
//   "parts of type piston"
//   "parts of type track"
//   "parts of type baron"                -- van berlo's wheel
//   "parts of type bonder"
//   "parts of type unbonder"
//   "parts of type bonder-prisma"        -- triplex bonder
//   "parts of type bonder-speed"         -- multi-bonder
//   "parts of type glyph-calcification"
//   "parts of type glyph-dispersion"
//   "parts of type glyph-disposal"
//   "parts of type glyph-duplication"
//   "parts of type glyph-life-and-death" -- glyph of animismus
//   "parts of type glyph-marker"         -- glyph of equilibrium
//   "parts of type glyph-projection"
//   "parts of type glyph-purification"
//   "parts of type glyph-unification"
//   "parts of type input"
//   "parts of type out-rep"              -- repeating (infinite) product
//   "parts of type out-std"              -- non-repeating product
//   "parts of type pipe"                 -- conduit
// the "number of track segments" metric iterates over each track piece, counts
//  all the segments in that piece, and returns the combined total.
// the "number of arms" metric returns the number of arms (the van berlo's wheel
//  is considered an arm).
// the "maximum absolute arm rotation" metric runs the solution to completion
//  while recording the absolute value of the rotation of each arm on every
//  cycle.  the maximum absolute rotation value is returned.
// the "overlap" metric places each part in the solution, counting each hex of
//  area that is placed on top of an existing part's area (not including
//  grabbers or arm bases that overlap track).
// the "duplicate reagents" metric counts how many reagents appear more than
//  once in the solution.  each appearance after the first of each reagent adds
//  one to the total count.  reagents that don't appear in the puzzle (and would
//  cause the solution to crash when loaded) are not counted.
// the "duplicate products" metric counts how many products appear more than
//  once in the solution in the same way that "duplicate reagents" counts
//  reagents.
// the "maximum track gap^2" metric iterates over each track piece, measuring
//  the distance between each segment in that piece and the next segment.  the
//  square of the maximum distance is returned.  this will always be either 1 or
//  0 for solutions created in-game.
// the "cost" metric returns the combined cost of all components.
// the "instructions" metric counts each instruction like the game does.
// the "instructions with hotkey <hotkeys>" metric counts instructions according
//  to their in-game hotkey (QWERTASDFG).  repeat and reset instructions are
//  expanded before instructions are counted -- specifying a C or V hotkey will
//  produce an error, since they have already been expanded.  you can specify
//  multiple hotkeys: e.g. "instructions with hotkey ADQE" to count rotations
//  and pivots.
// the "instruction tape period" metric counts how many cycles elapse before the
//  instruction tape loops.
// the "cycles" metric runs the solution to completion and returns the cycle the
//  last output was dropped.  if the solution doesn't run to completion, it
//  returns -1 -- run verifier_error() for the string describing the error.
// the "area (approximate)" metric runs the solution to completion in the same
//  way as "cycles", but returns the approximate area (to within about 1%).
// the "throughput cycles" and "throughput outputs" metrics run the solution
//  until it loops, then return the number of cycles in the loop and the number
//  of outputs produced.
// the "throughput cycles (unrestricted)" & "throughput outputs (unrestricted)"
//  metrics are similar to the above metrics, but give an approximate result for
//  some solutions that would be rejected.
// the "throughput waste" metric runs the solution until it loops.  if the
//  solution grows larger forever, a positive number is returned.  otherwise,
//  zero is returned.
// the "height" metric runs the solution to completion, then measures the length
//  of the area footprint perpendicular to each axis.  the smallest measured
//  length is returned.
// the "width*2" metric runs the solution to completion, then measures the
//  length of the area footprint parallel to each axis.  the smallest measured
//  length is multiplied by two (to avoid fractional results) and returned.
int verifier_evaluate_metric(void *verifier, const char *metric);

#endif
