from functools import total_ordering
from heapq import heapify, heappop
from random import random
from time import perf_counter_ns
import networkx

# graph is a list of lists


@total_ordering
class MergedVertex:
    base: int
    inner_vertices: list[int]
    phase_list_weight: float

    def __init__(self, vertices: list[int]):
        if len(vertices) == 0:
            raise ValueError("Supply at least one vertex!")
        self.base = vertices[0]
        self.inner_vertices = vertices
        self.phase_list_weight = -1

    def update_phase_list_weight(
        self, added_vertex: "MergedVertex", merged_weights: dict[int, dict[int, float]]
    ):
        if self.phase_list_weight < 0:
            self.phase_list_weight = merged_weights[added_vertex.base][self.base]
        else:
            self.phase_list_weight += merged_weights[added_vertex.base][self.base]

    def merge(self, other_vertex: "MergedVertex"):
        self.inner_vertices += other_vertex.inner_vertices

        return self

    def __lt__(self, other: "MergedVertex"):
        # we want the max!
        return self.phase_list_weight > other.phase_list_weight


def init_graph(n: int) -> list[MergedVertex]:
    return [MergedVertex([i]) for i in range(n)]


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
    cut_first_part: list[int] = second_to_last_vertex.inner_vertices.copy()
    for v in phase_list:
        cut_first_part.extend(v.inner_vertices)

    return cut_first_part


def create_phase_list(
    graph: list[MergedVertex],
    phase_list: list[MergedVertex],
    merged_weights: dict[int, dict[int, float]],
):
    while len(graph) > 2:
        # print([(m_v.phase_list_weight, m_v.base) for m_v in graph])
        next_vertex = heappop(graph)
        # print(next_vertex.base)
        phase_list.append(next_vertex)
        for v in graph:
            v.update_phase_list_weight(next_vertex, merged_weights)

    return graph, phase_list


def min_cut_phase(
    graph: list[MergedVertex],
    merged_weights: dict[int, dict[int, float]],
    active_bases: set[int],
):
    phase_list: list[MergedVertex] = []

    heapify(graph)

    graph, phase_list = create_phase_list(graph, phase_list, merged_weights)

    # now two vertices are left, which we will merge
    # first get second to last
    second_to_last_vertex = heappop(graph)
    # get the final one
    final_vertex = graph.pop()

    # update the final vertex, which will be the weight of the cut
    final_vertex.update_phase_list_weight(second_to_last_vertex, merged_weights)

    cut_first_part = create_first_cut(phase_list, second_to_last_vertex)

    cut = [cut_first_part, final_vertex.inner_vertices.copy()]
    cut_weight = final_vertex.phase_list_weight

    # merge the vertices
    final_vertex.merge(second_to_last_vertex)

    # uppdate the weights with the merged vertex
    merged_weights, active_bases = update_merged_weights(
        merged_weights, final_vertex.base, second_to_last_vertex.base, active_bases
    )

    return phase_list + [final_vertex], merged_weights, active_bases, cut, cut_weight


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

EdgeValue = tuple[float, list[tuple[int, int]]]


def sw_minimum_cut_phase(
    graph: list[list[EdgeValue]], a: int
) -> tuple[int, int, float, list[tuple[int, int]]]:
    graph_n = len(graph)
    A = [a]
    cut_edges: list[tuple[int, int]] = []
    max_cut_weight: float = -1

    while len(A) < graph_n:
        max_cut_weight = -1
        u = 0

        for v in range(graph_n):
            if v not in A:
                cut_weight = sum(graph[v][w][0] for w in A)
                if cut_weight > max_cut_weight:
                    max_cut_weight = cut_weight
                    u = v
                    cut_edges = []
                    for w in A:
                        cut_edges += graph[v][w][1]

        A.append(u)

    s = min(A[-1], A[-2])
    t = max(A[-1], A[-2])

    return s, t, max_cut_weight, cut_edges


def sw_minimum_cut(graph: list[list[EdgeValue]]):
    """Find the minimum cut of a graph using the Stoer-Wagner algorithm."""
    n = len(graph)

    min_cut = float("inf")
    contractions = []
    phase = 0
    best_edge_list = []

    while n > 1:
        a = 0  # Any vertex from V
        s, t, cut_weight, cut_edges = sw_minimum_cut_phase(graph[:n][:n], a)
        if cut_weight < min_cut:
            min_cut = cut_weight
            best_edge_list = cut_edges

        # Merge vertices s and t
        contractions.append((s, t))
        for i in range(n):
            if i != t:
                graph[s][i] = (
                    graph[s][i][0] + graph[t][i][0],
                    graph[s][i][1] + graph[t][i][1],
                )
                graph[i][s] = graph[s][i]

        for i in range(t, n - 1):
            for j in range(n):
                graph[i][j] = graph[i + 1][j]
                graph[j][i] = graph[i][j]

        n -= 1
        phase += 1

    return min_cut, best_edge_list


def edge_idx(lower_i: int, higher_j: int, n: int):
    """Returns the index of the edge using our edge-index convention."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)


def adjacency_matrix_from_vector(
    x_values: list[float], n: int
) -> tuple[list[list[EdgeValue]], list[list[float]]]:
    """Create an adjacency matrix from x_values"""
    if len(x_values) != (n * (n - 1)) // 2:
        raise ValueError("x_values does not have the correct length")

    adjacency_matrix = [[(0.0, [(i, j)]) for i in range(n)] for j in range(n)]
    weights: list[list[float]] = [[0 for j in range(n)] for i in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            idx = edge_idx(i, j, n)
            adjacency_matrix[i][j] = (x_values[idx], [(i, j)])
            adjacency_matrix[j][i] = adjacency_matrix[i][j]
            weights[i][j] = x_values[idx]
            weights[j][i] = x_values[idx]

    return adjacency_matrix, weights


def compute_min_cut(x_values: list[float], n: int):
    """Compute min cut using Stoer–Wagner algorithm. x_values should be 1D array with our edge
    index convention.
    """
    # we compute an adjacency matrix to more easily perform the stoer-wagner algorithm
    graph, weights = adjacency_matrix_from_vector(x_values, n)

    time_a = perf_counter_ns()

    # sw_min_cut, best_edge_list = sw_minimum_cut(graph)
    # we turn the list of edges back into array using our edge index convention
    # cut_edges = [edge_idx(e[0], e[1], n) for e in best_edge_list]

    time_b = perf_counter_ns()

    new_sw_min_cut, new_sw_min_cut_verts = min_cut(weights)

    time_c = perf_counter_ns()

    # print(sw_min_cut)
    print(new_sw_min_cut)

    # print(best_edge_list)
    # print(new_sw_min_cut_verts)

    print(
        f"Time first method {(time_b - time_a) / 1e6} ms. Second method {(time_c - time_b) / 1e6} ms"
    )

    # return min_cut, cut_edges


n = 200
m_edges = n * (n - 1) // 2

x_vector = [random() for _ in range(m_edges)]
# print(x_vector)

compute_min_cut(x_vector, n)
