from dataclasses import dataclass
import random
from typing import Any, Optional, Union

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


Edge = tuple[int, int]
VertTourNghbs = dict[int, tuple[int, int]]
VertNghbs = dict[int, set[int]]


def pick_x0(d_tour: VertTourNghbs, t0: int, costs: list[list[float]]):
    t0_nghbs = d_tour[t0]
    # first one is always t0
    return (t0, t0_nghbs[0]), (t0, t0_nghbs[1])


def exchange_gain(x_orig: Edge, y_repl: Edge, costs: list[list[float]]):
    orig_cost = costs[x_orig[0]][x_orig[1]]
    new_cost = costs[y_repl[0]][y_repl[1]]

    return orig_cost - new_cost


def pick_y0(n: int, d_tour: VertTourNghbs, t0: int, t1: int, costs: list[list[float]]) -> tuple[Optional[Edge], list[Edge], float]:
    in_tour = d_tour[t1]
    other_tour_nghb = in_tour[0] if in_tour[0] != t0 else in_tour[1]
    possible_y0 = [(t1, i) for i in range(n) if i != t0 and i != t1 and i != other_tour_nghb]

    x0 = (t0, t1)
    for maybe_y0 in possible_y0:
        gain = exchange_gain(x0, maybe_y0, costs)
        if gain > 0:
            possible_y0.remove(maybe_y0)
            return maybe_y0, possible_y0, gain

    return None, possible_y0, -1



# def tour_to_edge_matrix(n: int, tour: list[int]) -> list[list[float]]:
#     tour_edge_matrix = [[0 for _ in range(n)] for _ in range(n)]
#     for i in range(len(tour)):
#         tour_edge_matrix[tour[i]][tour[i+1]] = 1
#         tour_edge_matrix[tour[i+1]][tour[i]] = 1
    
#     tour_edge_matrix[tour[0]][tour[-1]] = 1
#     tour_edge_matrix[tour[-1]][tour[0]] = 1

#     return tour_edge_matrix


def d_tour_to_node_adj_dict(n: int, d_tour: VertTourNghbs) -> VertNghbs:
    node_edges = dict[int, set[int]]()
    
    for t in d_tour:
        node_edges[t] = set(d_tour[t])

    return node_edges


def is_node_edges_tour(n: int, node_edges: VertNghbs) -> tuple[bool, list[int]]:
    # maybe replace by whether t0 is reachable or something else
    tour = vert_ngbhs_to_tour(node_edges)
    return len(tour) == n, tour


def edge_eq(e1: Edge, e2: Edge) -> bool:
    return e1 == e2 or (e1[1], e1[0]) == e2


def edge_in_edge_list(e: Edge, edge_list: Union[list[Edge], set[Edge]]) -> bool:
    return e in edge_list or (e[1], e[0]) in edge_list


def copy_node_edges(node_edges: VertNghbs) -> VertNghbs:
    node_edge_copy = dict()
    for t in node_edges:
        node_edge_copy[t] = node_edges[t].copy()
    
    return node_edge_copy


def exchange_edges(node_edges: VertNghbs, x_remove: Edge, y_repl: Edge, copy=True) -> VertNghbs:
    if copy:
        node_edges = copy_node_edges(node_edges)
    t_x0, t_x1 = x_remove
    t_y0, t_y1 = y_repl


    node_edges[t_x0].remove(t_x1)
    node_edges[t_x1].remove(t_x0)

    node_edges[t_y0].add(t_y1)
    node_edges[t_y1].add(t_y0)

    return node_edges


def try_is_tour(n: int, tour_n_a_dict: VertNghbs, x_remove: Edge, y_repl: Edge) -> tuple[bool, list[int], VertNghbs]:
    maybe_tour_n_a_dict = exchange_edges(tour_n_a_dict, x_remove, y_repl)
    is_tour, maybe_tour = is_node_edges_tour(n, tour_n_a_dict)
    return is_tour, maybe_tour, maybe_tour_n_a_dict


def vert_ngbhs_to_tour(v_n: VertNghbs) -> list[int]:
    if len(v_n) == 0:
        return []
    current_node = next(iter(v_n))
    seen_nodes = {current_node}
    while True:
        adjacent = v_n[current_node]
        if len(adjacent) != 2:
            # not a valid tour
            return []
        not_seen = None
        for a in adjacent:
            if a not in seen_nodes:
                not_seen = a
                break
        if not_seen is None:
            break
        
        current_node = not_seen
        seen_nodes.add(not_seen)

    return list(seen_nodes)

@dataclass
class TourState:
    n: int
    t_base: int
    tour: list[int]
    d_tour: VertTourNghbs
    costs: list[list[float]]



def pick_xi(n: int, t_start: int, t_base: int, d_tour: VertTourNghbs, node_edges: VertNghbs, prev_ys: list[Edge]) -> Optional[tuple[Edge, Edge, list[int], VertNghbs]]:
    # node edges should have already removed x0, added y0, ... added yi-1
    t2i1_possible = d_tour[t_start]
    xi_0 = (t_start, t2i1_possible[0])
    xi_1 = (t_start, t2i1_possible[1])
    joined_0 = (t_base, t2i1_possible[0])
    joined_1 = (t_base, t2i1_possible[1])

    is_0_tour, maybe_tour, maybe_tour_n_a_dict = try_is_tour(n, node_edges, xi_0, joined_0)

    if not is_0_tour or edge_in_edge_list(xi_0, prev_ys):
        is_1_tour, maybe_tour, maybe_tour_n_a_dict = try_is_tour(n, node_edges, xi_1, joined_1)
        if not is_1_tour or edge_in_edge_list(xi_1, prev_ys):
            return None
        xi = xi_1
        joined = joined_1
        tour = maybe_tour
        tour_n_a_dict = maybe_tour_n_a_dict
    else:
        xi = xi_0
        joined = joined_0
        tour = maybe_tour
        tour_n_a_dict = maybe_tour_n_a_dict

    # TODO is this guaranteed?
    # if xi in prev_ys:
    #     return None
    
    return xi, joined, tour, tour_n_a_dict


def pick_yi_xi1(n: int, t2i1: int, t0: int, xi: Edge, partial_gain: float, d_tour: VertTourNghbs, tour_n_a_dict: VertNghbs, costs: list[list[float]], prev_xs: list[Edge], prev_ys: list[Edge]) -> Optional[tuple[Edge, Edge, Edge, float, list[int], VertNghbs]]:
    tour_nghbs = d_tour[t2i1]

    possible_yi = [(t2i1, i) for i in range(n) if i != t2i1 and i not in tour_nghbs]

    for maybe_yi in possible_yi:

        gain_i = exchange_gain(xi, maybe_yi, costs)
        if partial_gain + gain_i > 0:
            maybe_yi_rev = (maybe_yi[1], maybe_yi[0])
            if edge_in_edge_list(maybe_yi, prev_xs) or edge_in_edge_list(maybe_yi_rev, prev_xs):
                continue
            t2i2 = maybe_yi[1]
            updated_tour_n_a_dict = exchange_edges(tour_n_a_dict, xi, maybe_yi)
            maybe_prev_ys = prev_ys.copy() + [maybe_yi]
            maybe_xi1 = pick_xi(n, t2i2, t0, d_tour, updated_tour_n_a_dict, maybe_prev_ys)
            if maybe_xi1 is None:
                continue
            xi1, joined_i1, new_tour, updated_tour_n_a_dict = maybe_xi1
            
            possible_yi.remove(maybe_yi)
            return maybe_yi, xi1, joined_i1, gain_i, new_tour, updated_tour_n_a_dict

    return None
    
    

# def tour_to_edges(tour: list[int]) -> list[tuple[int, int]]:
#     edges = []
#     for i in range(len(tour)):
#         edges.append((tour[i], tour[i+1]))
#     edges.append((tour[-1], tour[0]))

#     return edges

def next_t2i(i: int, t_chain: list[int]):
    t_chain.append(0)      
    t2i = t_chain[2*i]

    return t2i, t_chain

UsedEdges = dict[int, set[Edge]]
SelectedEdges = dict[int, Edge]


def next_x(i: int, used_x: UsedEdges, sel_y: SelectedEdges, cur_graph: VertNghbs, ts: TourState):
    if i == 0:
        t_start = ts.t_base
    else:
        yim1 = sel_y[i-1]
        t_start = yim1[1]

    base_nghbs = ts.d_tour[t_start]
    for nghb in base_nghbs:
        maybe_x = (t_start, nghb)
        if edge_in_edge_list(maybe_x, used_x.get(0, [])):
            continue
        
        if i >= 1:
            joined = (ts.t_base, nghb)
            is_tour, maybe_tour, updated_tgraph = try_is_tour(ts.n, cur_graph, maybe_x, joined)

            if not is_tour:
                continue

            if edge_in_edge_list(maybe_x, list(sel_y.values())):
                continue
        else:
            # no y_i to replace with
            updated_tgraph = cur_graph
            maybe_tour = ts.tour

        return maybe_x, maybe_tour, updated_tgraph

    return None
        
    
    


def next_y(i: int, used_y: UsedEdges):
    pass


def add_create_edge(i: int, dct: dict[int, set[Edge]], e: Edge):
    if i in dct:
        dct[i].add(e)
    else:
        dct[i] = {e}


def evaluate_x(i: int, used_x: UsedEdges, used_y: UsedEdges, sel_x: SelectedEdges, sel_y: SelectedEdges, cur_graph: VertNghbs, ts: TourState):
    x_i_result = next_x(i, used_x, sel_y, cur_graph, ts)
    if x_i_result is not None:
        x_i, new_tour, updated_tgraph = x_i_result
        # do check if better
        # else
        sel_x[i] = x_i
        add_create_edge(i, used_x, x_i)
        return evaluate_y(i, used_x, used_y, sel_x, sel_y)
    
    if i == 1:
        used_x.pop(1, None)
        used_y.pop(0, None)
    elif i == 0:
        # exhausted all x_0 so go to different t_base
        return None
    else:
        raise ValueError("x_i should exist for higher i!")
    



def evaluate_y(i: int, used_x: UsedEdges, used_y: UsedEdges, sel_x: SelectedEdges, sel_y: SelectedEdges):
    y_i_result = next_y(i, used_y)
    y_i = (0, 0) # TODO change
    if y_i_result is not None:   
        sel_y[i] = y_i
        add_create_edge(i, used_y, y_i)
        i = i+1
        return evaluate_x(i, used_x, used_y, sel_x, sel_y)
    if i >= 2:
        # since we just tried for y_1 we know no more exist
        i = 1
        return evaluate_x(1, used_x, used_y, sel_x, sel_y)
    elif i == 1:
        used_y.pop(1, None)
        return evaluate_x(1, used_x, used_y, sel_x, sel_y)
    elif i == 0:
        used_y.pop(0, None)
        return evaluate_x(0, used_x, used_y, sel_x, sel_y)


def try_t0(n: int, t0: int, tour: list[int], costs: list[list[float]]) -> list[int]:
    d_tour = dict_tour(tour)
    tour_n_a_dict = d_tour_to_node_adj_dict(n, d_tour)
    i = 0
    t_chain = [t0]
    while True:
        t_chain.append(0)
        x0, unpicked_x0 = pick_x0(d_tour, t0, costs)
        prev_xs = [x0]
        # first one is t0
        t1 = x0[1]
        t_chain[1] = t1
        y0, unpicked_y0, gain_0 = pick_y0(n, d_tour, t0, t1, costs)
        partial_gain_0 = gain_0
        
        if y0 is None:
            # try another t0 (break from inner while loop)
            break

        t_chain[2] = y0[1]
        prev_ys = [y0]
        updated_tour_n_a_dict = exchange_edges(tour_n_a_dict, x0, y0)
        i = 1
        t2i, t_chain = next_t2i(i, t_chain)
        maybe_x1 = pick_xi(n, t2i, t0, d_tour, updated_tour_n_a_dict, prev_ys)
        if maybe_x1 is None:
            # TODO relax the tour having to be valid
            raise ValueError("x1 one should exist!")
        
        x1, joined_1, new_tour, updated_tour_n_a_dict = maybe_x1
        x1_joined_gain = exchange_gain(x1, joined_1, costs)
        if partial_gain_0 + x1_joined_gain > 0:
            return new_tour
        prev_xs.append(x1)
        xi = x1
        partial_gain_im1 = partial_gain_0
        while True:
            t2i1 = t_chain[2*i+1]
            maybe_yi_xi1 = pick_yi_xi1(n, t2i1, t0, xi, partial_gain_im1, d_tour, tour_n_a_dict, costs, prev_xs, prev_ys)
            
            if maybe_yi_xi1 is None:
                break

            yi, xi1, joined_i1, gain_i, new_tour, updated_tour_n_a_dict = maybe_yi_xi1
            partial_gain_i = partial_gain_im1 + gain_i
            t_chain += [0, 0]
            t_chain[2*i+2] = yi[1]
            t_chain[2*i+3] = xi1

            i = i + 1
            partial_gain_im1 = partial_gain_i
            
            
            # t2i = t_chain[2*i]
            # maybe_yi = pick_xi(n, t2i, t0, d_tour, updated_tour_n_a_dict, prev_ys)
            
            

            # t2i1 = xi[1]
            # t_chain[2*i + 1] = t2i1

def lin_kernighan(n: int, costs: list[list[float]]):
    tour = random_tour(n)

    while True:
        untried_t0 = tour.copy()

        better_tour = None

        while len(untried_t0) > 0:
            t0 = untried_t0.pop()
            tour_t0, is_better = try_t0(n, t0, tour, costs)

            if is_better:
                better_tour = tour_t0

        if better_tour is None:
            break
    
    
    

    






