import random
from typing import Optional

def dict_tour(tour: list[int]) -> dict[int, tuple[int, int]]:
    # adjacent edges in the tour
    d: dict[int, tuple[int, int]] = dict()
    final_i = len(tour)-1
    for i in range(len(tour)):
        if i == 0:
            other_i1 = 1
            other_i2 = final_i
        elif i == final_i:
            other_i1 = 0
            other_i2 = final_i-1
        else:
            other_i1 = i+1
            other_i2 = i-1
        
        d[tour[i]] = (tour[other_i1], tour[other_i2])
    
    return d

def random_tour(n: int) -> list[int]:
    randoms = [(random.random(), i) for i in range(n)]
    randoms.sort(key=lambda l: l[0])
    tour = list(map(lambda r: r[1], randoms))
    return tour


def pick_x0(d_tour: dict[int, tuple[int, int]], t0: int, costs: list[list[int]]):
    t0_nghbs = d_tour[t0]
    # first one is always t0
    return (t0, t0_nghbs[0]), (t0, t0_nghbs[1])


def exchange_gain(x_orig: tuple[int, int], y_repl: tuple[int, int], costs: list[list[int]]):
    orig_cost = costs[x_orig[0]][x_orig[1]]
    new_cost = costs[y_repl[0]][y_repl[1]]

    return orig_cost - new_cost


def pick_y0(n: int, d_tour: dict[int, tuple[int, int]], t0: int, t1: int, costs: list[list[int]]) -> tuple[Optional[tuple[int, int]], list[tuple[int, int]]]:
    in_tour = d_tour[t1]
    other_tour_nghb = in_tour[0] if in_tour[0] != t0 else in_tour[1]
    possible_y0 = [(t1, i) for i in range(n) if i != t0 and i != t1 and i != other_tour_nghb]

    x0 = (t0, t1)
    for maybe_y0 in possible_y0:
        gain = exchange_gain(x0, maybe_y0, costs)
        if gain > 0:
            possible_y0.remove(maybe_y0)
            return maybe_y0, possible_y0

    return None, possible_y0



def tour_to_edge_matrix(n: int, tour: list[int]) -> list[list[int]]:
    tour_edge_matrix = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(len(tour)):
        tour_edge_matrix[tour[i]][tour[i+1]] = 1
        tour_edge_matrix[tour[i+1]][tour[i]] = 1
    
    tour_edge_matrix[tour[0]][tour[-1]] = 1
    tour_edge_matrix[tour[-1]][tour[0]] = 1

    return tour_edge_matrix


def d_tour_to_node_edges(n: int, d_tour: dict[int, tuple[int, int]]) -> dict[int, set[int]]:
    node_edges = dict[int, set[int]]()
    
    for t in d_tour:
        node_edges[t] = set(d_tour[t])

    return node_edges


def is_node_edges_tour(n: int, node_edges: dict[int, set[int]]) -> bool:
    # maybe replace by whether t0 is reachable or something else
    if len(node_edges) == 0:
        return False
    current_node = next(iter(node_edges))
    seen_nodes = {current_node}
    while True:
        adjacent = node_edges[current_node]
        if len(adjacent) != 2:
            return False
        not_seen = None
        for a in adjacent:
            if a not in seen_nodes:
                not_seen = a
                break
        if not_seen is None:
            break
        
        current_node = not_seen
        seen_nodes.add(not_seen)

    return len(seen_nodes) == n


def edge_eq(e1: tuple[int, int], e2: tuple[int, int]):
    return e1 == e2 or (e1[1], e1[0]) == e2


def copy_node_edges(node_edges: dict[int, set[int]]):
    node_edge_copy = dict()
    for t in node_edges:
        node_edge_copy[t] = node_edges[t].copy()
    
    return node_edge_copy

def pick_xi(n: int, t2i: int, t0: int, d_tour: dict[int, tuple[int, int]], node_edges: dict[int, set[int]], prev_ys: list[tuple[int, int]]):
    # node edges should have already removed x0, added y0, ... added yi-1
    t2i1_possible = d_tour[t2i]

    for t2i1 in t2i1_possible:
        n_edges = copy_node_edges(node_edges)
    
        # remove x_i
        n_edges[t2i].remove(t2i1)
        n_edges[t2i1].remove(t2i)

        # add (t0, t_{2i+1})
        n_edges[t0].add(t2i1)
        n_edges[t2i1].add(t0)

    is_tour = is_node_edges_tour(n, node_edges)


    

def tour_to_edges(tour: list[int]) -> list[tuple[int, int]]:
    edges = []
    for i in range(len(tour)):
        edges.append((tour[i], tour[i+1]))
    edges.append((tour[-1], tour[0]))

    return edges


def lin_kernighan(n: int, costs: list[list[int]]):
    tour = random_tour(n)
    
    
    untried_t0 = tour.copy()

    while True:
        d_tour = dict_tour(tour)
        i = 0
        t0 = untried_t0.pop()
        while True:
            x0, unpicked_x0 = pick_x0(d_tour, t0, costs)
            # first one is t0
            t1 = x0[1]
            y0, unpicked_y0 = pick_y0(n, d_tour, t0, t1, costs)

            if y0 is None:
                # try another t0 (break from inner while loop)
                break
            
            i = i + 1




