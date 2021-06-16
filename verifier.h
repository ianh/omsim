#ifndef OM_VERIFIER_H
#define OM_VERIFIER_H

// this "verifier" API is designed to be called from an FFI.

// takes the path of a .puzzle and a .solution file.  returns a pointer to a
// newly-allocated verifier object.  if the files couldn't be parsed, the
// verifier_error() function will return a non-null string value describing the
// error.
void *verifier_create(const char *puzzle_filename, const char *solution_filename);

// destroys a verifier object created by verifier_create(), freeing its memory.
void verifier_destroy(void *verifier);

// returns a string describing any error that occurred during solution parsing,
// decoding, or verification.  if there haven't been any errors, returns a null
// pointer.  returned strings are valid for the lifetime of the program and
// don't need to be freed.
const char *verifier_error(void *verifier);

// verifier_evaluate_metric() evaluates a metric for the given solution.
// the "parsed cycles", "parsed cost", "parsed area", and "parsed instructions"
//  metrics return the value recorded in the solution file.  if no value is
//  recorded, returns zero.
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
// the "cost" metric returns the combined cost of all components.
// the "instructions" metric counts each instruction like the game does.
// the "cycles" metric runs the solution to completion and returns the cycle the
//  last output was dropped.  if the solution doesn't run to completion, it
//  returns -1 -- run verifier_error() for the string describing the error.
// the "area (approximate)" metric runs the solution to completion in the same
//  way as "cycles", but returns the approximate area (to within about 1%).
// the "throughput cycles" and "throughput outputs" metrics run the solution
//  until it loops, then return the number of cycles in the loop and the number
//  of outputs produced.
// the "height" metric runs the solution to completion, then measures the length
//  of the area footprint perpendicular to each axis.  the smallest measured
//  length is returned.
// the "width*2" metric runs the solution to completion, then measures the
//  length of the area footprint parallel to each axis.  the smallest measured
//  length is multiplied by two (to avoid fractional results) and returned.
int verifier_evaluate_metric(void *verifier, const char *metric);

#endif
