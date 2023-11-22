import numpy as np

graph = [
    [0, 12, 0, 0, 0, 0, 12],
    [12, 0, 8, 12, 0, 0, 0],
    [0, 8, 0, 11, 3, 0, 9],
    [0, 12, 11, 0, 11, 10, 13],
    [0, 0, 3, 11, 0, 6, 7],
    [0, 0, 0, 10, 6, 0, 9],
    [12, 0, 9, 13, 7, 9, 0],
]

n = len(graph)
m_edges = (n * (n - 1)) // 2


def length_tour(graph, tour):
    # kan sneller
    length = 0  # graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        length += graph[tour[i]][tour[i + 1]]

    return length


tour = [0, 1, 2, 3, 4, 5, 6, 0]
print(length_tour(graph, tour))


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


def simulated_annealing(graph, T, r):
    return


print(tour)
