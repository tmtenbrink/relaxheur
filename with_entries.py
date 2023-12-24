from functools import total_ordering
from heapq import heapify, heappop, heappush
from random import random
from time import perf_counter_ns
import networkx

from fib_heap import FibonacciHeap
from heap_dct import heapdict

# graph is a list of lists


class MergedVertex:
    base: int
    inner_vertices: list[int]

    def __init__(self, vertices: list[int], base=None):
        if len(vertices) == 0:
            raise ValueError("Supply at least one vertex!")

        self.base = vertices[0]
        self.inner_vertices = vertices

    def merge(self, other_vertex: "MergedVertex"):
        self.inner_vertices += other_vertex.inner_vertices

        return self


def new_phase_light_weight(
    prev_value: float,
    base: int,
    added_base: int,
    merged_weights: dict[int, dict[int, float]],
):
    if prev_value < 0:
        return merged_weights[added_base][base]
    else:
        return prev_value + merged_weights[added_base][base]


@total_ordering
class HeapVertexEntry:
    v: MergedVertex
    valid: bool
    phase_list_weight: float

    def __init__(self, v: MergedVertex, phase_light_weight: float):
        self.v = v
        self.phase_list_weight = phase_light_weight
        self.valid = True

    def updated_entry(
        self, added_base: int, merged_weights: dict[int, dict[int, float]]
    ):
        updated_phase_light_weight = new_phase_light_weight(
            self.phase_list_weight, self.v.base, added_base, merged_weights
        )
        self.valid = False
        return HeapVertexEntry(self.v, updated_phase_light_weight)

    def __lt__(self, other: "HeapVertexEntry"):
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
    heap: list[HeapVertexEntry],
    phase_list: list[MergedVertex],
    merged_weights: dict[int, dict[int, float]],
):
    print(heap)
    for _ in range(len(heap) - 2):
        # print_heap_values(heap)
        # print([(m_v.phase_list_weight, m_v.base) for m_v in graph])
        next_vertex, _ = pop_from_heap(heap)
        # print(next_vertex.base)
        # print(f"{[(m_v.phase_list_weight, m_v.base) for m_v in graph]} after pop")
        # print(next_vertex.base)
        phase_list.append(next_vertex)
        update_phase_list_weight(heap, next_vertex, merged_weights)

    return heap, phase_list


def print_heap_values(heap):
    _, m_v_dict = heap

    print([(v[0], v[1][1].key) for v in m_v_dict.items()])


def update_phase_list_weight(
    heap: list[HeapVertexEntry],
    next_vertex: MergedVertex,
    merged_weights: dict[int, dict[int, float]],
):
    for h_e in heap:
        new_entry = h_e.updated_entry(next_vertex.base, merged_weights)
        heappush(heap, new_entry)

    # fib_heap, m_v_dict = heap

    # for v_and_nd in m_v_dict.values():
    #     v, nd = v_and_nd
    #     new_weight = v.update_phase_list_weight(next_vertex, merged_weights)
    #     fib_heap.decrease_key(nd, -new_weight)

    # heap_dct, m_v_dict = heap
    # for v in m_v_dict.values():
    #     new_weight = v.update_phase_list_weight(next_vertex, merged_weights)
    #     heap_dct[v.base] = -new_weight

    # heapify(heap)


# def heap_len(heap) -> int:
#     return len(heap)

#     # _, m_v_dict = heap

#     # return len(m_v_dict)


def create_heap(graph: list[MergedVertex]):
    return [HeapVertexEntry(v, -1) for v in graph]

    # fib_heap = FibonacciHeap()
    # m_v_dict = {}
    # for v in graph:
    #     fib_nd = fib_heap.insert(-v.phase_list_weight, v.base)
    #     m_v_dict[v.base] = (v, fib_nd)

    # return (fib_heap, m_v_dict)

    # heap = heapdict()
    # m_v_dict = {}
    # for v in graph:
    #     heap[v.base] = -v.phase_list_weight
    #     m_v_dict[v.base] = v

    # return (heap, m_v_dict)


def pop_from_heap(heap: list[HeapVertexEntry]) -> tuple[MergedVertex, float]:
    entry = heappop(heap)
    while not entry.valid:
        entry = heappop(heap)

    return entry.v, entry.phase_list_weight

    # fib_heap, m_v_dict = heap

    # v_base = fib_heap.extract_min().value

    # return m_v_dict.pop(v_base)[0]

    # heap_dct, m_v_dict = heap

    # v_base, _ = heap_dct.popitem()

    # return m_v_dict.pop(v_base)


def min_cut_phase(
    graph: list[MergedVertex],
    merged_weights: dict[int, dict[int, float]],
    active_bases: set[int],
):
    phase_list: list[MergedVertex] = []

    graph_heap = create_heap(graph)

    graph_heap, phase_list = create_phase_list(graph_heap, phase_list, merged_weights)

    # now two vertices are left, which we will merge
    # first get second to last
    second_to_last_vertex, _ = pop_from_heap(graph_heap)
    # get the final one
    final_vertex, final_prev_weight = pop_from_heap(graph_heap)

    # update the final vertex, which will be the weight of the cut
    cut_weight = new_phase_light_weight(
        final_prev_weight, final_vertex.base, second_to_last_vertex.base, merged_weights
    )

    cut_first_part = create_first_cut(phase_list, second_to_last_vertex)

    cut = [cut_first_part, final_vertex.inner_vertices.copy()]

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


n = 3
m_edges = n * (n - 1) // 2

x_vector = [random() for _ in range(m_edges)]
# print(x_vector)

compute_min_cut(x_vector, n)
