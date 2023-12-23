

from typing import Optional


EdgeValue = tuple[float, list[tuple[int, int]]]
Edge = tuple[int, int]
VertTourNghbs = dict[int, tuple[int, int]]
# all neighbors of each node
VertNghbs = dict[int, set[int]]

UsedEdges = dict[int, set[Edge]]
SelectedEdges = dict[int, Edge]


def try_is_tour(
    n: int, tour_n_a_dict: VertNghbs, x_remove: Edge, y_repl: Edge, moves: list[tuple[Edge, Edge]], base_tgraph: VertNghbs
) -> tuple[bool, list[int]]:

    is_tour, maybe_tour = is_node_edges_tour(n, tour_n_a_dict, x_remove, y_repl)

    return is_tour, maybe_tour


def exchange_edges(
    node_edges: VertNghbs,
    x_remove: Edge,
    y_repl: Edge,
    copy=True,
) -> VertNghbs:
    if copy:
        node_edges = copy_node_edges(node_edges)
    t_x0, t_x1 = x_remove
    t_y0, t_y1 = y_repl

    edges_x0 = node_edges[t_x0]
    edges_x1 = node_edges[t_x1]
    edges_x0.discard(t_x1)
    edges_x1.discard(t_x0)
    node_edges[t_y0].add(t_y1)
    node_edges[t_y1].add(t_y0)

    return node_edges


def vert_ngbhs_to_tour_seen(v_n: VertNghbs) -> list[int]:
    if len(v_n) == 0:
        return []
    current_node = next(iter(v_n))
    first_node = current_node
    prev_node = None
    tour = [current_node]
    while True:
        adjacent = v_n[current_node]
        if len(adjacent) != 2:
            # not a valid tour
            return []
        not_seen = None
        for adj in adjacent:
            if adj != prev_node:
                not_seen = adj
                break
        if not_seen is None:
            return tour
        
        if not_seen == first_node:
            return tour

        prev_node = current_node
        current_node = not_seen
        
        tour.append(not_seen)


def is_node_edges_tour(n: int, node_edges: VertNghbs, x_remove: Edge, y_repl: Edge) -> tuple[bool, list[int]]:
    maybe_tour_n_a_dict = exchange_edges(
        node_edges, x_remove, y_repl, copy=False
    )
    tour = vert_ngbhs_to_tour_seen(maybe_tour_n_a_dict)
    exchange_edges(
        maybe_tour_n_a_dict, y_repl, x_remove, copy=False
    )

    if len(tour) != n:
        return False, []
    return True, tour


def copy_node_edges(node_edges: VertNghbs) -> VertNghbs:
    node_edge_copy = dict()
    for t in node_edges:
        node_edge_copy[t] = node_edges[t].copy()

    return node_edge_copy