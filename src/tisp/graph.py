from tisp.types import BestNeighbors, Costs, EdgesByIndex, VertexEdges


def edge_idx(lower_i: int, higher_j: int, n: int):
    """Returns the index of the edge using our edge-index convention."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)


def get_edge_idxs_for_v(vertex: int, n: int):
    """Return the indexes of all edges that contain the vertex."""
    lower_edges = [edge_idx(other_v, vertex, n) for other_v in range(0, vertex)]
    higher_edges = [edge_idx(vertex, other_v, n) for other_v in range(vertex + 1, n)]

    return lower_edges + higher_edges


def edge_idxs_for_all_v(n: int) -> VertexEdges:
    return [get_edge_idxs_for_v(v, n) for v in range(n)]


def get_edges_by_index(n: int) -> EdgesByIndex:
    edges_by_index: dict[int, tuple[int, int]] = {}
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            edges_by_index[index] = (i, j)

            index += 1

    return edges_by_index


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


def copy_costs(costs: Costs):
    return [row.copy() for row in costs]


def best_neighbors(costs: list[list[float]], group_len=5) -> BestNeighbors:
    n = len(costs)
    best_nbs: list[list[int]] = []
    for i, cost_row in enumerate(costs):
        best_costs = sorted(zip(cost_row, range(n)))
        best_group = []
        for _, j in best_costs:
            if i == j:
                continue
            best_group.append(j)
            if len(best_group) == group_len:
                break
        best_nbs.append(best_group)
    return best_nbs
