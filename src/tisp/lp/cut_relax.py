from numbers import Real
from typing import cast
import mip
from tisp.graph import adjacency_matrix_from_vector, edge_idx

from tisp.types import MergedVertex, EdgeValues, LPConstants, LPModel


EPSILON = 0.0000001


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

    while len(graph) > 2:
        next_vertex = graph.pop(next_vertex_idx)
        phase_list.append(next_vertex)

        added_vertex_weights = merged_weights[next_vertex[1]]

        highest_wght = 0
        next_vertex_idx = 0

        i = 0
        for prev_value, base, inner_verts in graph:

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

    second_to_last_vertex = graph.pop(penultimate_idx)

    final_vertex = graph.pop()

    penultimate_base = second_to_last_vertex[1]

    cut_first_part = create_first_cut(phase_list, second_to_last_vertex)

    cut = [cut_first_part, final_vertex[2].copy()]
    cut_weight = new_phase_list_weight(
        final_vertex[0], final_vertex[1], merged_weights[penultimate_base]
    )

    final_vertex = merge_vertex(final_vertex, second_to_last_vertex)

    merged_weights, active_bases = update_merged_weights(
        merged_weights, final_vertex[1], penultimate_base, active_bases
    )

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

    min_cut_value = [[], []]
    min_cut_weight = float("inf")

    while len(graph) > 1:
        graph, merged_weights, active_bases, cut, cut_weight = min_cut_phase(
            graph, merged_weights, active_bases
        )

        if cut_weight < min_cut_weight:
            min_cut_value = cut
            min_cut_weight = cut_weight

    return min_cut_weight, min_cut_value


def compute_min_cut(x_values: list[float], n: int) -> tuple[float, list[int]]:
    """Compute min cut using Stoer–Wagner algorithm. x_values should be 1D array with our edge
    index convention.
    """
    # we compute an adjacency matrix to more easily perform the stoer-wagner algorithm
    graph = adjacency_matrix_from_vector(x_values, n)
    min_cut_weight, cut = min_cut(graph)

    cut_edges = []
    for n1 in cut[0]:
        for n2 in cut[1]:
            edge = edge_idx(n1, n2, n) if n1 < n2 else edge_idx(n2, n1, n)
            cut_edges.append(edge)

    return min_cut_weight, cut_edges


def separation(
    model: mip.Model, edge_values: EdgeValues, edge_vars: list[mip.Var], c: LPConstants
) -> bool:
    """Tests whether x is in P_subtour, if not it adds the violated inequality to the model.
    1) The constraint >=0 is already defined in the char_vector.
    2) We check if x(delta(v))=2, for all v, with some tolerance epsilon.
    3) We check x(delta(U))>= 2, for all U, by finding the minimum cut and checking if it
    is larger than 2 (with tolerance epsilon). Note that if the minimimum cut is larger
    than 2, we know that all cuts are larger than 2.

    It is easy to see that constraints 1 and 2 are checked in polynomial time. Constraint 3
    has exponentially many inequalities (as there are exponentially many U), but can be checked
    in polynomial time since the min-cut can be found in polynomial time by the Stoer-Wagner algorithm.

    Therefore, our separation algorithm is polynomial time.
    """

    n, _, vertex_edges, _, _ = c

    for v in range(n):
        # the columns are the vertices
        edges = vertex_edges[v]
        x_sum = 0
        for e in edges:
            x_sum += edge_values[e]

        if abs(x_sum - 2) > EPSILON:
            constr_vars = dict([(edge_vars[e], cast(Real, 1.0)) for e in edges])
            model += mip.LinExpr(expr=constr_vars) == 2, f"x(delta({v}))==2"
            return False

    cut_weight, min_cut_edges = compute_min_cut(edge_values, n)
    if cut_weight < 2 - EPSILON:
        subtour_vars = dict([(edge_vars[e], cast(Real, 1.0)) for e in min_cut_edges])
        model += mip.LinExpr(expr=subtour_vars) >= cast(Real, 2.0), f"x(delta(cut)>=2"

        return False

    return True
