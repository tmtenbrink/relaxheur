

EdgeValue = tuple[float, list[tuple[int, int]]]
Edge = tuple[int, int]
VertTourNghbs = dict[int, tuple[int, int]]
VertNghbs = dict[int, set[int]]

UsedEdges = dict[int, set[Edge]]
SelectedEdges = dict[int, Edge]


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


def try_is_tour(
    n: int, tour_n_a_dict: VertNghbs, x_remove: Edge, y_repl: Edge
) -> tuple[bool, list[int]]:
    maybe_tour_n_a_dict = exchange_edges(
        tour_n_a_dict, x_remove, y_repl, require_existence=False, copy=False
    )
    is_tour, maybe_tour = is_node_edges_tour(n, maybe_tour_n_a_dict)
    maybe_tour_n_a_dict = exchange_edges(
        tour_n_a_dict, y_repl, x_remove, require_existence=False, copy=False
    )
    return is_tour, maybe_tour


# def try_is_tour(
#     n: int, tour_n_a_dict: VertNghbs, x_remove: Edge, y_repl: Edge
# ) -> tuple[bool, list[int]]:
#     maybe_tour_n_a_dict = exchange_edges(
#         tour_n_a_dict, x_remove, y_repl, require_existence=False, copy=True
#     )
#     is_tour, maybe_tour = is_node_edges_tour(n, maybe_tour_n_a_dict)
#     # exchange_edges(
#     #     maybe_tour_n_a_dict, y_repl, x_remove, require_existence=False, copy=False
#     # )
#     return is_tour, maybe_tour



def exchange_edges(
    node_edges: VertNghbs,
    x_remove: Edge,
    y_repl: Edge,
    copy=True,
    require_existence=True,
) -> VertNghbs:
    if copy:
        node_edges = copy_node_edges(node_edges)
    t_x0, t_x1 = x_remove
    t_y0, t_y1 = y_repl

    edges_x0 = node_edges[t_x0]
    edges_x1 = node_edges[t_x1]
    if require_existence or t_x1 in edges_x0:
        edges_x0.remove(t_x1)
    if require_existence or t_x0 in edges_x1:
        edges_x1.remove(t_x0)

    node_edges[t_y0].add(t_y1)
    node_edges[t_y1].add(t_y0)

    return node_edges


def vert_ngbhs_to_tour(v_n: VertNghbs) -> list[int]:
    if len(v_n) == 0:
        return []
    current_node = next(iter(v_n))
    seen_nodes = {current_node}
    seen_order = [current_node]
    while True:
        adjacent = v_n[current_node]
        if len(adjacent) != 2:
            # not a valid tour
            return []
        not_seen = None
        for adj in adjacent:
            if adj not in seen_nodes:
                not_seen = adj
                break
        if not_seen is None:
            break

        current_node = not_seen
        seen_nodes.add(not_seen)
        seen_order.append(not_seen)

    return seen_order


def is_node_edges_tour(n: int, node_edges: VertNghbs) -> tuple[bool, list[int]]:
    # maybe replace by whether t0 is reachable or something else
    tour = vert_ngbhs_to_tour(node_edges)
    # return len(tour) == n, tour
    if len(tour) != n:
        return False, []
    return True, tour


def copy_node_edges(node_edges: VertNghbs) -> VertNghbs:
    node_edge_copy = dict()
    for t in node_edges:
        node_edge_copy[t] = node_edges[t].copy()

    return node_edge_copy