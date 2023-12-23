from numbers import Real
from typing import cast
import mip
from tisp.graph import adjacency_matrix_from_vector, edge_idx

from tisp.types import EdgeValue, EdgeValues, LPConstants, LPModel


EPSILON = 0.0000001


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


def compute_min_cut(x_values: list[float], n: int) -> tuple[float, list[int]]:
    """Compute min cut using Stoer–Wagner algorithm. x_values should be 1D array with our edge
    index convention.
    """
    # we compute an adjacency matrix to more easily perform the stoer-wagner algorithm
    graph = adjacency_matrix_from_vector(x_values, n)
    min_cut, best_edge_list = sw_minimum_cut(graph)
    # we turn the list of edges back into array using our edge index convention
    cut_edges = [edge_idx(e[0], e[1], n) for e in best_edge_list]

    return min_cut, cut_edges



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
    