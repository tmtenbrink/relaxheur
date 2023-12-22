

from typing import Optional


EdgeValue = tuple[float, list[tuple[int, int]]]
Edge = tuple[int, int]
VertTourNghbs = dict[int, tuple[int, int]]
# all neighbors of each node
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


def maybe_tour_is_tour(n: int, tour_n_a_dict: VertNghbs, x_remove: Edge, y_repl: Edge):
    maybe_tour = vert_nghbs_to_tour_exchange(n, tour_n_a_dict, x_remove, y_repl)
    if len(maybe_tour) != n:
        return False, []
    return True, maybe_tour


def try_is_tour(
    n: int, tour_n_a_dict: VertNghbs, x_remove: Edge, y_repl: Edge
) -> tuple[bool, list[int]]:
    is_tour, maybe_tour = maybe_tour_is_tour(n, tour_n_a_dict, x_remove, y_repl)
    # if len(tour_new_mthd) != n:
    #     return False, []
    # if len(tour_new_mthd) != n:
    #     return False, []
    # return True, tour
    
    # is_tour, maybe_tour = is_node_edges_tour(n, tour_n_a_dict, x_remove, y_repl)
    # if len(maybe_tour) != len(tour_new_mthd):
    #     tour_new_mthd = vert_nghbs_to_tour_exchange(tour_n_a_dict, x_remove, y_repl)
    #     print("Bad!")
    # we exchange back to not modify the original list
    
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


# def remove_other(node: int, removed_e: Edge, added_e: Edge, neighbors: set[int], seen_nodes: set[int]) -> Optional[int]:
#     to_skip = []
#     to_add = []
#     if node == removed_e[0]:
#         to_skip = removed_e[1]
#     elif node == added_e[1]:
#         to_skip = removed_e[1]

#     if node == added_e:


def other_node(node: int, e: Edge):
    if node == e[0]:
        return e[1]
    elif node == e[1]:
        return e[0]
    else:
        return None


def vert_nghbs_to_tour_exchange(n: int, v_n: VertNghbs, x_remove: Edge, y_repl: Edge):
    if len(v_n) == 0:
        return []
    
    current_node = next(iter(v_n))
    first_node = current_node
    
    prev_node = None
    tour = [current_node]

    while True:
        remov_adj = other_node(current_node, x_remove)
        added_adj = other_node(current_node, y_repl)
        added_adj = [added_adj] if added_adj is not None and added_adj != prev_node and added_adj != remov_adj else []
        unseen_nghbs = [nb for nb in v_n[current_node] if nb != prev_node and nb != remov_adj] + added_adj

        # if we are the first node there should be two possibilities
        required_length = 1 if prev_node is not None else 2

        if len(unseen_nghbs) == 0:
            return tour

        if len(unseen_nghbs) != required_length:
            return []
        
        next_nb = unseen_nghbs[0]

        if next_nb == first_node:
            return tour

        tour.append(next_nb)
        prev_node = current_node
        current_node = next_nb


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


def is_node_edges_tour(n: int, node_edges: VertNghbs, x_remove: Edge, y_repl: Edge) -> tuple[bool, list[int]]:
    maybe_tour_n_a_dict = exchange_edges(
        node_edges, x_remove, y_repl, require_existence=False, copy=False
    )
    
    # maybe replace by whether t0 is reachable or something else
    tour = vert_ngbhs_to_tour(maybe_tour_n_a_dict)

    exchange_edges(
        maybe_tour_n_a_dict, y_repl, x_remove, require_existence=False, copy=False
    )

    if len(tour) != n:
        return False, []
    return True, tour


def copy_node_edges(node_edges: VertNghbs) -> VertNghbs:
    node_edge_copy = dict()
    for t in node_edges:
        node_edge_copy[t] = node_edges[t].copy()

    return node_edge_copy