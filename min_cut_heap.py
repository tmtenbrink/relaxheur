from functools import total_ordering
from heapq import heapify, heappop
from random import random
from time import perf_counter_ns
import networkx

# graph is a list of lists


MergedVertex = tuple[float, int, list[int]]


def merge_vertex(base: MergedVertex, other: MergedVertex) -> MergedVertex:
    return 0, base[1], base[2] + other[2]


def new_phase_list_weight(
    prev_value: float, base: int, added_vertex_weights: dict[int, float]
):
    return prev_value + added_vertex_weights[base]


def init_graph(n: int) -> list[MergedVertex]:
    return [(0, i, [i]) for i in range(n)]


def update_merged_weights(
    merged_weights: dict[int, dict[int, float]],
    vertex_base: int,
    other_vertex_base: int,
    active_bases: set[int],
):
    active_bases.difference_update([vertex_base, other_vertex_base])

    weights = merged_weights[other_vertex_base]
    other_weights = merged_weights[vertex_base]
    merged_weights[vertex_base] = {
        i: weights[i] + other_weights[i] for i in active_bases
    }
    for v in active_bases:
        merged_weights[v][vertex_base] = weights[v] + other_weights[v]

    active_bases.add(vertex_base)

    return merged_weights, active_bases


def create_first_cut(
    phase_list: list[MergedVertex], second_to_last_vertex: MergedVertex
):
    # we haven't added this to the phase list as we don't want to include it in the new graph
    # the cut will be from the final vertex to the rest, so we do want to add the vertices
    cut_first_part: list[int] = second_to_last_vertex[2].copy()
    for v in phase_list:
        cut_first_part.extend(v[2])

    return cut_first_part


def create_phase_list(
    graph: list[MergedVertex],
    phase_list: list[MergedVertex],
    merged_weights: dict[int, dict[int, float]],
):
    next_vertex_idx = 0

    # print(f"starting phase with graph {graph}")
    # print(f"merged weights are {merged_weights}")

    while len(graph) > 2:
        next_vertex = graph.pop(next_vertex_idx)
        phase_list.append(next_vertex)
        # print(f"next vertex is {next_vertex}")
        # print(f"remaining graph is {graph}")

        added_vertex_weights = merged_weights[next_vertex[1]]
        # print(f"weights of next vertex is {added_vertex_weights}")

        highest_wght = 0
        next_vertex_idx = 0
        # # print(graph)
        i = 0
        for prev_value, base, inner_verts in graph:
            # print(f"updating remaining vertex {base}")
            ph_lst_weight = new_phase_list_weight(
                prev_value, base, added_vertex_weights
            )
            graph[i] = (ph_lst_weight, base, inner_verts)

            if ph_lst_weight > highest_wght:
                next_vertex_idx = i
                highest_wght = ph_lst_weight

            i += 1

    return graph, phase_list, next_vertex_idx


def reset_phase_weights(graph: list[MergedVertex]) -> list[MergedVertex]:
    return [(0, v[1], v[2]) for v in graph]


def min_cut_phase(
    graph: list[MergedVertex],
    merged_weights: dict[int, dict[int, float]],
    active_bases: set[int],
):
    phase_list: list[MergedVertex] = []

    graph, phase_list, penultimate_idx = create_phase_list(
        graph, phase_list, merged_weights
    )

    # print(f"remaining graph after phase list: {graph}")
    # now two vertices are left, which we will merge
    # first get second to last
    second_to_last_vertex = graph.pop(penultimate_idx)
    # get the final one
    final_vertex = graph.pop()
    # SOMETHING WRONG WITH WHICH ONES ARE MERGED
    # print(f"second_to_last: {second_to_last_vertex}; final: {final_vertex}")

    # update the final vertex, which will be the weight of the cut
    penultimate_base = second_to_last_vertex[1]

    cut_first_part = create_first_cut(phase_list, second_to_last_vertex)

    cut = [cut_first_part, final_vertex[2].copy()]
    cut_weight = new_phase_list_weight(
        final_vertex[0], final_vertex[1], merged_weights[penultimate_base]
    )

    # merge the vertices
    final_vertex = merge_vertex(final_vertex, second_to_last_vertex)

    # print(f"merging vertices {penultimate_base} and {final_vertex[1]}")

    # uppdate the weights with the merged vertex
    # print(f"weights before merging: {merged_weights}")
    # print(f"active bases: {active_bases}")
    merged_weights, active_bases = update_merged_weights(
        merged_weights, final_vertex[1], penultimate_base, active_bases
    )
    # print(f"active bases after update: {active_bases}")

    return (
        reset_phase_weights(phase_list) + [final_vertex],
        merged_weights,
        active_bases,
        cut,
        cut_weight,
    )


def min_cut(weights: list[list[float]]):
    n = len(weights)
    graph = init_graph(n)
    merged_weights = {i: {j: weights[i][j] for j in range(n)} for i in range(n)}
    active_bases = set(range(n))

    min_cut_value = None
    min_cut_weight = float("inf")

    while len(graph) > 1:
        graph, merged_weights, active_bases, cut, cut_weight = min_cut_phase(
            graph, merged_weights, active_bases
        )

        if cut_weight < min_cut_weight:
            min_cut_value = cut
            min_cut_weight = cut_weight

    return min_cut_weight, min_cut_value


EPSILON = 0.0000001


def edge_idx(lower_i: int, higher_j: int, n: int):
    """Returns the index of the edge using our edge-index convention."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)


def adjacency_matrix_from_vector(x_values: list[float], n: int) -> list[list[float]]:
    """Create an adjacency matrix from x_values"""
    if len(x_values) != (n * (n - 1)) // 2:
        raise ValueError("x_values does not have the correct length")

    weights: list[list[float]] = [[0 for j in range(n)] for i in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            idx = edge_idx(i, j, n)
            weights[i][j] = x_values[idx]
            weights[j][i] = x_values[idx]

    return weights


def compute_min_cut(x_values: list[float], n: int):
    """Compute min cut using Stoer–Wagner algorithm. x_values should be 1D array with our edge
    index convention.
    """
    # we compute an adjacency matrix to more easily perform the stoer-wagner algorithm
    weights = adjacency_matrix_from_vector(x_values, n)

    # sw_min_cut, best_edge_list = sw_minimum_cut(graph)
    # we turn the list of edges back into array using our edge index convention
    # cut_edges = [edge_idx(e[0], e[1], n) for e in best_edge_list]

    time_b = perf_counter_ns()

    new_sw_min_cut, new_sw_min_cut_verts = min_cut(weights)

    time_c = perf_counter_ns()

    # # print(sw_min_cut)
    print(new_sw_min_cut)

    # # print(best_edge_list)
    print(new_sw_min_cut_verts)

    print(f"Time {(time_c - time_b) / 1e6} ms")

    # return min_cut, cut_edges


n = 200
m_edges = n * (n - 1) // 2

x_vector = [random() * (0 if random() < 0.999 else 1) for _ in range(m_edges)]
print(x_vector)
# # print(x_vector)

compute_min_cut(x_vector, n)
