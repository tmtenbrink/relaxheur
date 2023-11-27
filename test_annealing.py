import argparse
from pathlib import Path
import random
import math
import numpy as np


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
    if tour == None:
        return 0

    length = graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        length += graph[tour[i]][tour[i + 1]]

    return length


def nearest_neighbor1(graph, starting_node=None):
    # gets error if distances are not unique
    if starting_node == None:
        starting_node = random.randrange(0, len(graph), 1)
    tour = [starting_node]
    for i in range(len(graph) - 1):
        row = graph[tour[i]]
        filtered_values = [d for d in row if d > 0 and row.index(d) not in tour]
        min_distance = min(filtered_values)
        next_node = row.index(min_distance)
        tour.append(next_node)

    tour.append(tour[starting_node])

    return tour


def nearest_neighbor2(graph, starting_node=None):
    # fixed non-uniqueness of distances
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


def nearest_neighbor3(graph, max_retries=10):
    for _ in range(max_retries):
        starting_node = random.randrange(0, len(graph), 1)
        tour = [starting_node]  # Choose random starting node
        for i in range(len(graph) - 1):
            row = graph[tour[i]]
            filtered_values = [d for d in row if d > 0 and row.index(d) not in tour]

            if not filtered_values:
                # If filtered_values is empty, start over with a new random starting node
                break
            else:
                min_distance = min(filtered_values)
                next_node = row.index(min_distance)

            tour.append(next_node)

        tour.append(tour[starting_node])

        # If the tour length is equal to the number of nodes, a valid tour is found
        if len(tour) == len(graph) + 1:
            return tour

    # If no valid tour is found after the maximum number of retries, return None
    print("No tour found")
    return None


def simulated_annealing(graph, T, r, L=10000, max_no_improvement=500):
    S = nearest_neighbor1(graph)  # set initial solution S to be greedy solution
    epsilon = 1e-6
    no_improvement_count = 0

    while T > epsilon:  # while not frozen
        for _ in range(L):
            if no_improvement_count >= max_no_improvement:
                no_improvement_count = 0
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

    n, graph_l = parse_as_adj_matrix(inst_path)

    for i in range(n):
        greedy_tour2 = nearest_neighbor2(graph_l, i)
        print(i, length_tour(graph_l, greedy_tour2))
        try:
            greedy_tour1 = nearest_neighbor1(graph_l, i)
            print(length_tour(graph_l, greedy_tour1))
        except:
            print("fail")

    # annealing_tour = simulated_annealing(graph_l, T=100, r=0.9)
    # print(annealing_tour)
    # print(length_tour(graph_l, annealing_tour))


if __name__ == "__main__":
    run()
