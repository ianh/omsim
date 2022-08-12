#ifndef OM_VERIFIER_H
#define OM_VERIFIER_H

// this "verifier" API is designed to be called from an FFI.

// find the puzzle name within the bytes of a solution.  returns a pointer to
// the first byte of the puzzle name within the solution_bytes, or a null
// pointer if the solution failed to parse.  if successful, also sets the value
// pointed to by name_length to the length of the puzzle name.
const char *verifier_find_puzzle_name_in_solution_bytes(const char *solution_bytes,
 int solution_length, int *name_length);

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

// set whether the verifier fails if an output of the correct shape but with the
// wrong atoms is dropped on the output with puzzle file index output_index.
void verifier_set_fails_on_wrong_output(void *verifier, int output_index, int fails_on_wrong_output);
// get the puzzle file index of the current wrong output or -1 if no wrong
// outputs were detected.
int verifier_wrong_output_index(void *verifier);
// get the atom type at the given offset in the current wrong output.  the u and
// v offsets are as they appear in the puzzle file.  if there's no wrong output
// or no atom at this position, the value -1 is returned.
int verifier_wrong_output_atom(void *verifier, int u, int v);
// clears the current wrong output so that verifier_wrong_output_index will
// return -1 until another wrong output is detected.
void verifier_wrong_output_clear(void *verifier);

// returns a string describing any error that occurred during solution parsing,
// decoding, or verification.  if there haven't been any errors, returns a null
// pointer.  returned strings are valid for the lifetime of the program and
// don't need to be freed.  calling this function doesn't clear the error;
// to clear the error, create a new verifier or call verifier_error_clear().
const char *verifier_error(void *verifier);

// these functions return more information about errors occurring during the
// simulation.
int verifier_error_cycle(void *verifier);
int verifier_error_location_u(void *verifier);
int verifier_error_location_v(void *verifier);

// after calling this function, verifier_error() will return null again.
void verifier_error_clear(void *verifier);

// verifier_evaluate_metric() evaluates a metric for the given solution.  if
// there is an error, verifier_evaluate_metric() returns -1 and sets an error
// string, which verifier_error() will return.
//
// for an up-to date list of metrics, see the "simulator metrics" section of
// <http://events.critelli.technology/static/metrics.html>.
int verifier_evaluate_metric(void *verifier, const char *metric);

// used to compute the lexicographic cycles metric.  the first interval is
// the cycle on which the first output is dropped; the second interval is the
// difference between the cycles on which the second and the first inputs are
// dropped; the third interval is the difference between the third and the
// second; and so on through verifier_number_of_output_intervals().
// if the index returned by verifier_output_intervals_repeat_after() is
// non-negative, the output intervals from that index onward repeat cyclically.
// the output is specified using the puzzle file's output index.
int verifier_number_of_output_intervals(void *verifier, int which_output);
int verifier_output_interval(void *verifier, int which_output, int which_interval);
int verifier_output_intervals_repeat_after(void *verifier, int which_output);

#endif
