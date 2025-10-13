#include "parse.h"
#include "verifier.h"

#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_usage(FILE *f)
{
    fprintf(f, "Usage: omsim [--metric <metric name>]* [--disable-limits] "
               "[--puzzle-file <file> | --puzzle-folder <folder>] solution(s)...\n");
}

void print_help(FILE *f)
{
    print_usage(f);
    fprintf(f, "Options:\n");
    fprintf(f, "  -p --puzzle-file <file>     Specify a puzzle file\n");
    fprintf(f, "  -f --puzzle-folder <folder> Specify a folder containing puzzles\n");
    fprintf(f, "  -m --metric <name>          Evaluate the given metric\n");
    fprintf(f, "  -s --steady-state           Simulate up to steady-state\n");
    fprintf(f, "  -d --disable-limits         Disable simulation limits\n");
    fprintf(f, "  -h --help                   Show this help message and exit\n");
}

void print_verifier_error(void *verifier)
{
    if (verifier_error(verifier)) {
        fputs(verifier_error(verifier), stdout);
        if (verifier_error_cycle(verifier) != 0)
            printf(" on cycle %d at %d %d", verifier_error_cycle(verifier),
                   verifier_error_location_u(verifier), verifier_error_location_v(verifier));
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    const char *puzzle_file = NULL;
    const char *puzzle_folder = NULL;

    const char *requested_metrics[64];
    unsigned int number_of_requested_metrics = 0;

    bool steady_state = false;
    bool disable_limits = false;

    char *const *solution_paths = NULL;
    unsigned int number_of_solution_paths = 0;

    static struct option long_options[] = {{"puzzle-file", required_argument, NULL, 'p'},
                                           {"puzzle-folder", required_argument, NULL, 'f'},
                                           {"metric", required_argument, NULL, 'm'},
                                           {"steady-state", no_argument, NULL, 's'},
                                           {"disable-limits", no_argument, NULL, 'd'},
                                           {"help", no_argument, NULL, 'h'},
                                           {NULL, 0, NULL, 0}};

    int opt;
    while ((opt = getopt_long(argc, argv, "p:f:m:sdh", long_options, NULL)) != -1) {
        switch (opt) {
        case 'm':
            if (number_of_requested_metrics >=
                (int)(sizeof(requested_metrics) / sizeof(requested_metrics[0]))) {
                fprintf(stderr, "too many metrics\n");
                return -1;
            }
            requested_metrics[number_of_requested_metrics++] = optarg;
            break;
        case 'd':
            disable_limits = true;
            break;
        case 'p':
            puzzle_file = optarg;
            break;
        case 'f':
            puzzle_folder = optarg;
            break;
        case 's':
            steady_state = true;
            break;
        case 'h':
            print_help(stdout);
            return 0;
        case '?':
        default:
            print_help(stderr);
            return -1;
        }
    }

    // If puzzle_file or puzzle_folder is set, use them instead of positional puzzle_path
    if (puzzle_file && puzzle_folder) {
        fprintf(stderr, "cannot specify both -p|--puzzle-file and -f|--puzzle-folder\n");
        return -1;
    } else if (!puzzle_file && !puzzle_folder) {
        fprintf(stderr, "must specify either -p|--puzzle-file or -f|--puzzle-folder\n");
        return -1;
    }

    number_of_solution_paths = argc - optind;
    if (number_of_solution_paths < 1) {
        print_usage(stderr);
        return -1;
    }
    solution_paths = argv + optind;

    int retval = 0;
    for (int i = 0; i < number_of_solution_paths; ++i) {
        const char *solution_path = solution_paths[i];
        char *puzzle_path = NULL; // must be freed
        if (puzzle_file) {
            puzzle_path = strdup(puzzle_file);
        } else if (puzzle_folder) {
            struct solution_file *sf = parse_solution_file(solution_path);
            if (sf) {
                size_t puzzle_path_length =
                    strlen(puzzle_folder) + 1 + sf->puzzle.length + strlen(".puzzle") + 1;
                puzzle_path = malloc(puzzle_path_length);
                if (!puzzle_path) {
                    fprintf(stderr, "out of memory\n");
                    return -1;
                }
                snprintf(puzzle_path, puzzle_path_length, "%s/%.*s.puzzle", puzzle_folder,
                         (int)sf->puzzle.length, sf->puzzle.bytes);
                free_solution_file(sf);
            } else {
                fprintf(stderr, "couldn't find puzzle name in solution file %s\n", solution_path);
                retval = -1;
                continue;
            }
        }

        void *verifier = verifier_create(puzzle_path, solution_path);
        if (verifier_error(verifier)) {
            print_verifier_error(verifier);
            retval = -1;
            goto cleanup;
        }
        if (disable_limits)
            verifier_disable_limits(verifier);

        if (number_of_solution_paths > 1)
            printf("%s:\n", solution_path);

        if (number_of_requested_metrics > 0) {
            for (int i = 0; i < number_of_requested_metrics; ++i) {
                double value = verifier_evaluate_approximate_metric(verifier, requested_metrics[i]);
                printf("%s: ", requested_metrics[i]);
                if (verifier_error(verifier)) {
                    print_verifier_error(verifier);
                    retval = -1;
                    verifier_error_clear(verifier);
                } else
                    printf("%g\n", value);
            }
        } else {
            int cost = verifier_evaluate_metric(verifier, "cost");
            if (verifier_error(verifier)) {
                print_verifier_error(verifier);
                retval = -1;
                goto cleanup;
            }
            int instructions = verifier_evaluate_metric(verifier, "instructions");
            printf("%dg/%di@0", cost, instructions);

            int cycles = verifier_evaluate_metric(verifier, "cycles");
            if (verifier_error(verifier)) {
                printf("\n");
                print_verifier_error(verifier);
                retval = -1;
                goto cleanup;
            }
            int area = verifier_evaluate_metric(verifier, "area");
            printf(" %dc/%da@V", cycles, area);

            if (steady_state) {
                int outINF = verifier_evaluate_metric(verifier, "per repetition outputs");
                if (verifier_error(verifier)) {
                    printf("\n");
                    print_verifier_error(verifier);
                } else {
                    double rate =
                        (double)verifier_evaluate_metric(verifier, "per repetition cycles") /
                        outINF;
                    printf(" %gr@INF\n", rate);
                }
            } else
                printf("\n");
        }

    cleanup:
        verifier_destroy(verifier);
        free(puzzle_path);
    }
    return retval;
}
