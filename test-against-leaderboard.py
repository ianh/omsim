from ctypes import *
import json
from math import ceil
from pathlib import Path
import signal
import sys

if len(sys.argv) < 3:
    print('usage: test-against-leaderboard.py [puzzle directory] [leaderboard directory]')
    exit(1)

def measure_rate(m):
    cycles = m("throughput cycles")
    outputs = m("throughput outputs")
    if outputs <= 0 or cycles < 0:
        return None
    else:
        return ceil(100 * cycles / outputs) / 100

def check_rate_metric(m, val):
    if val < 0 or not measure_rate(m):
        return None
    else:
        return val

leaderboard_metrics = {
    "cost": lambda m: m("cost"),
    "instructions": lambda m: m("instructions"),
    "overlap": lambda m: m("overlap") != 0,
    "trackless": lambda m: m("number of track segments") == 0,
    "cycles": lambda m: m("cycles"),
    "area": lambda m: m("area"),
    "width": lambda m: m("width*2") / 2,
    "height": lambda m: m("height"),
    "rate": measure_rate,
    "areaINF": lambda m: check_rate_metric(m, m("steady state area")),
    "heightINF": lambda m: check_rate_metric(m, m("steady state height")),
    "widthINF": lambda m: check_rate_metric(m, m("steady state width*2") / 2),
}

signal.signal(signal.SIGINT, signal.SIG_DFL)
lv = cdll.LoadLibrary("libverify.so")
lv.verifier_create.restype = c_void_p
lv.verifier_error.restype = c_char_p
lv.verifier_find_puzzle_name_in_solution_bytes.restype = POINTER(c_char)

puzzles = {}
for path in Path(sys.argv[1]).rglob('*.puzzle'):
    puzzles[path.stem] = path

leaderboard = Path(sys.argv[2])
for path in leaderboard.rglob('*.json'):
    with path.open() as f:
        metadata = json.load(f)
    solution_path = leaderboard / metadata["dataPath"]
    solution_bytes = solution_path.read_bytes()
    puzzle_name_length = c_int()
    puzzle_name_string = lv.verifier_find_puzzle_name_in_solution_bytes(c_char_p(solution_bytes), c_int(len(solution_bytes)), byref(puzzle_name_length))
    puzzle_path = puzzles[puzzle_name_string[:puzzle_name_length.value].decode('utf-8')]
    v = lv.verifier_create(c_char_p(bytes(puzzle_path)), c_char_p(bytes(solution_path)))
    for metric, expected in metadata["score"].items():
        # if the leaderboard doesn't track width/height for a puzzle, don't validate it.
        if (metric == "height" or metric == "heightINF") and not metadata["score"]["height"]:
            continue
        if (metric == "width" or metric == "widthINF") and not metadata["score"]["width"]:
            continue
        f = leaderboard_metrics[metric]
        measured = f(lambda sim_metric: lv.verifier_evaluate_metric(c_void_p(v), c_char_p(sim_metric.encode('utf-8'))))
        int_expected = None
        try:
            int_expected = int(expected)
        except:
            pass
        if measured != expected and measured != int_expected:
            print(puzzle_path, solution_path, metric, measured, expected)
    lv.verifier_destroy(c_void_p(v))
