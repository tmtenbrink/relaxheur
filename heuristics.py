import argparse
from pathlib import Path
import random
import math
from time import perf_counter_ns
from typing import Literal, Optional, Union, overload

from linkernighan import lin_kernighan


def parse_line(ln: str, as_float=False) -> Union[list[float], list[int]]:
    convert = float if as_float else int
    return list(map(lambda i: convert(i), ln.rstrip().split(" ")))


def get_inst_path():
    parser = argparse.ArgumentParser(description="Process some integers.")

    parser.add_argument("inst_name", help="Filename of instance")

    args = parser.parse_args()

    return Path(args.inst_name)


@overload
def parse_as_adj_matrix(inst_path: Path, as_float: Literal[True]) -> list[list[float]]:
    ...


@overload
def parse_as_adj_matrix(inst_path: Path, as_float: Literal[False]) -> list[list[int]]:
    ...


@overload
def parse_as_adj_matrix(
    inst_path: Path,
) -> list[list[int]]:
    ...


def parse_as_adj_matrix(inst_path: Path, as_float=False):  # type: ignore
    with open(inst_path, "r") as f:
        lines = f.readlines()

    # line_0 = parse_line(lines[0])
    adj_matrix = list(map(lambda ln: parse_line(ln, as_float), lines[1:]))
    # n = line_0[0]

    return adj_matrix


def length_tour(graph, tour):
    # kan sneller
    length = graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        length += graph[tour[i]][tour[i + 1]]

    return length


def nearest_neighbor(graph, starting_node=None):
    if starting_node == None:
        starting_node = random.randrange(0, len(graph), 1)
    tour = [starting_node]
    for i in range(len(graph) - 1):
        row = graph[tour[i]]
        row_new = [0 if i in tour else val for i, val in enumerate(row)]
        filtered_values = [d for d in row_new if d > 0]
        min_distance = min(filtered_values)
        next_node = row_new.index(min_distance)
        tour.append(next_node)

    tour.append(tour[starting_node])

    return tour


def simulated_annealing(graph, T, r, L=1000, max_no_improvement=1e9):
    S = nearest_neighbor(graph)  # set initial solution S to be greedy solution
    epsilon = 1e-6

    while T > epsilon:  # while not frozen
        no_improvement_count = 0
        for _ in range(L):
            if no_improvement_count >= max_no_improvement:
                break
            # Pick random S_new in neighborhood of S:
            # random start and end point sub-tour reversal
            S_new = S.copy()
            i = random.choice(range(1, len(S) - 2))
            j = random.choice(range(i + 1, len(S) - 1))
            subtour = S_new[i : j + 1]
            S_new[i : j + 1] = subtour[::-1]

            # Check if S_new is better than S or not
            diff = length_tour(graph, S_new) - length_tour(graph, S)
            if diff <= 0:
                S = S_new
                no_improvement_count = 0
            else:
                prob = math.exp(-diff / T)
                S = random.choices((S, S_new), weights=(1 - prob, prob))[0]
                no_improvement_count += 1
        T *= r
    return S


def run():
    # inst_path = get_inst_path()
    # inst_path = Path("tsp/burma14.dat")
    inst_path = Path("tsp/gr48.dat")

    graph_l = parse_as_adj_matrix(inst_path)

    print("Value of heuristic")

    # # Greedy Heuristic: Nearest neighbor
    start_time = perf_counter_ns()
    greedy_tour = nearest_neighbor(graph_l)
    greedy_time = (perf_counter_ns() - start_time) / 1e9
    greedy_length = length_tour(graph_l, greedy_tour)
    print(f"- Nearest Neighbor: {greedy_length}, {greedy_time}s")

    # # Simulated Annealing1
    start_time = perf_counter_ns()
    annealing_tour = simulated_annealing(graph_l, T=100, r=0.95, L=1000)
    annealing_time = (perf_counter_ns() - start_time) / 1e9
    annealing_length = length_tour(graph_l, annealing_tour)
    print(f"- Simulated Annealing: {annealing_length}, {annealing_time}s  (L=1000)")

    # Simulated Annealing2
    start_time = perf_counter_ns()
    annealing_tour = simulated_annealing(
        graph_l, T=100, r=0.95, L=10000, max_no_improvement=100
    )
    annealing_time = (perf_counter_ns() - start_time) / 1e9
    annealing_length = length_tour(graph_l, annealing_tour)
    print(
        f"- Simulated Annealing: {annealing_length}, {annealing_time}s  (L=10000, max_no_improvement=100)"
    )

    # Lin-Kernighan
    graph_l_fl = parse_as_adj_matrix(inst_path, as_float=True)
    start_time = perf_counter_ns()
    lin_tour, gain = lin_kernighan(graph_l_fl)
    lin_time = (perf_counter_ns() - start_time) / 1e9
    lin_length = length_tour(graph_l, lin_tour)
    print(f"- Lin-Kernighan: {lin_length}, {lin_time}s")


if __name__ == "__main__":
    run()
