import argparse
from pathlib import Path
import random
import math


def parse_line(ln: str):
    return list(map(lambda i: int(i), ln.rstrip().split(" ")))


def get_inst_path():
    parser = argparse.ArgumentParser(description="Process some integers.")

    parser.add_argument("inst_name", help="Filename of instance")

    args = parser.parse_args()

    return Path(args.inst_name)


def parse_as_adj_matrix(inst_path: Path):
    with open(inst_path, "r") as f:
        lines = f.readlines()

    line_0 = parse_line(lines[0])
    n = line_0[0]
    adj_matrix = list(map(lambda ln: parse_line(ln), lines[1:]))

    return n, adj_matrix


def length_tour(graph, tour):
    # kan sneller
    length = graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        length += graph[tour[i]][tour[i + 1]]

    return length


def nearest_neighbor(graph):
    tour = [0]
    for i in range(len(graph) - 1):
        row = graph[tour[i]]
        filtered_values = [d for d in row if d > 0 and row.index(d) not in tour]
        min_distance = min(filtered_values)
        next_node = row.index(min_distance)
        tour.append(next_node)

    tour.append(tour[0])

    return tour


def simulated_annealing(graph, T, r, L=100):
    S = nearest_neighbor(graph)  # set initial solution S to be greedy solution
    epsilon = 1e-6

    while T > epsilon:  # while not frozen
        for _ in range(L):
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
            else:
                prob = math.exp(-diff / T)
                S = random.choices((S, S_new), weights=(1 - prob, prob))[0]
        T *= r
    return S


def run():
    # inst_path = get_inst_path()
    # inst_path = Path("tsp/burma14.dat")
    inst_path = Path("tsp/att48.dat")

    n, graph_l = parse_as_adj_matrix(inst_path)

    greedy_tour = nearest_neighbor(graph_l)
    print(greedy_tour)
    print(length_tour(graph_l, greedy_tour))

    annealing_tour = simulated_annealing(graph_l, T=100, r=0.9, L=1000)
    print(annealing_tour)
    print(length_tour(graph_l, annealing_tour))


if __name__ == "__main__":
    run()
