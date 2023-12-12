import argparse
from pathlib import Path
import random
import math
from time import perf_counter_ns
from functools import reduce
from dataclasses import dataclass
import random
from typing import Iterable, Optional, Union, Literal, Optional, Union, overload


@overload
def parse_line(ln: str, as_float: Literal[True]) -> list[float]:
    ...


@overload
def parse_line(ln: str) -> list[int]:
    ...



def parse_line(ln: str, as_float=False) -> Union[list[float], list[int]]:
    convert = float if as_float else int
    return list(map(lambda i: convert(i), ln.rstrip().split(" ")))


def get_inst_path():
    parser = argparse.ArgumentParser(description="Process some integers.")

    parser.add_argument("inst_name", help="Filename of instance")

    args = parser.parse_args()

    return Path(args.inst_name)


@overload
def parse_as_adj_matrix(inst_path: Path, as_float: Literal[True]) -> list[list[float]]:
    ...


@overload
def parse_as_adj_matrix(
    inst_path: Path,
) -> list[list[int]]:
    ...


def parse_as_adj_matrix(inst_path: Path, as_float=False):
    with open(inst_path, "r") as f:
        lines = f.readlines()

    # line_0 = parse_line(lines[0])
    if as_float:
        adj_matrix = list(map(lambda ln: parse_line(ln, True), lines[1:]))
    else:
        adj_matrix = list(map(lambda ln: parse_line(ln), lines[1:]))
    # n = line_0[0]

    return adj_matrix


def length_tour(graph, tour):
    length = graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        length += graph[tour[i]][tour[i + 1]]

    return length



def nearest_neighbor(graph: list[list[int]], starting_node=None):
    if starting_node is None:
        starting_node = random.randrange(0, len(graph), 1)

    n = len(graph)
    remaining_nodes = list(range(n))
    cur_node = starting_node
    remaining_nodes.remove(starting_node)
    tour = [cur_node]
    # we need n-1 more nodes after starting with the random node
    for _ in range(n-1):
        # get the costs for all nodes
        costs_row = graph[cur_node]
        # get a start index
        cur_node = remaining_nodes[0]
        for n in remaining_nodes[1:]:
            if costs_row[n] < costs_row[cur_node]:
                cur_node = n
        tour.append(cur_node)
        remaining_nodes.remove(cur_node)

    return tour


def simulated_annealing(graph, T, r, L=1000, max_no_improvement=1e9):
    S = nearest_neighbor(graph)  # set initial solution S to be greedy solution
    epsilon = 1e-6

    while T > epsilon:  # while not frozen
        no_improvement_count = 0
        for _ in range(L):
            if no_improvement_count >= max_no_improvement:
                break
            # Pick random S_new in neighborhood of S:
            # random start and end point sub-tour reversal
            S_new = S.copy()
            i = random.choice(range(1, len(S) - 2))
            j = random.choice(range(i + 1, len(S) - 1))
            subtour = S_new[i : j + 1]
            S_new[i : j + 1] = subtour[::-1]

            # Check if S_new is better than S or not
            diff = length_tour(graph, S_new) - length_tour(graph, S)
            if diff <= 0:
                S = S_new
                no_improvement_count = 0
            else:
                prob = math.exp(-diff / T)
                S = random.choices((S, S_new), weights=(1 - prob, prob))[0]
                no_improvement_count += 1
        T *= r
    return S




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
    return shuffle_iter(range(n))


def shuffle_iter(itrble: Iterable[int]) -> list[int]:
    """Returns a shuffled copy of the iterable as a list."""
    randoms = [(random.random(), i) for i in itrble]
    randoms.sort(key=lambda l: l[0])
    shuffled = list(map(lambda r: r[1], randoms))
    return shuffled


Edge = tuple[int, int]
VertTourNghbs = dict[int, tuple[int, int]]
VertNghbs = dict[int, set[int]]

UsedEdges = dict[int, set[Edge]]
SelectedEdges = dict[int, Edge]


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


def edge_in_edge_list(e: Edge, edge_list: Iterable[Edge]) -> bool:
    return e in edge_list or (e[1], e[0]) in edge_list


def copy_node_edges(node_edges: VertNghbs) -> VertNghbs:
    node_edge_copy = dict()
    for t in node_edges:
        node_edge_copy[t] = node_edges[t].copy()
    
    return node_edge_copy


def exchange_edges(node_edges: VertNghbs, x_remove: Edge, y_repl: Edge, copy=True, require_existence=True) -> VertNghbs:
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

def try_is_tour(n: int, tour_n_a_dict: VertNghbs, x_remove: Edge, y_repl: Edge) -> tuple[bool, list[int]]:
    maybe_tour_n_a_dict = exchange_edges(tour_n_a_dict, x_remove, y_repl, require_existence=False)
    is_tour, maybe_tour = is_node_edges_tour(n, maybe_tour_n_a_dict)
    return is_tour, maybe_tour


@dataclass
class BaseTour:
    n: int
    tbase: int
    tour: list[int]
    d_tour: VertTourNghbs
    tgraph: VertNghbs
    costs: list[list[float]]


XResult = Optional[tuple[Edge, list[int], Edge]]
Tours = dict[int, tuple[VertNghbs, float]]
TourResult = tuple[list[int], float]

# Note that this object is used contextually (i.e. whether it is for next_x, next_y, etc.)
@dataclass
class Iteration:
    i: int
    next_sel: Literal['x', 'y']
    # Whether in previous step it was already generated (only applicable to xi)
    already_xi: Optional[XResult]
    used_x: UsedEdges
    used_y: UsedEdges
    sel_x: SelectedEdges
    sel_y: SelectedEdges
    tours: Tours
    bt: BaseTour
    result: Optional[TourResult]
    
    def with_i(self, i: int):
        self.i = i
        return self
    
    def with_nxt(self, nxt: Literal['x', 'y']):
        self.next_sel = nxt
        return self

    def var_next_x(self):
        return self.i, self.bt, self.sel_y, self.used_x, self.tours
    
    def var_next_y(self):
        return self.i, self.bt, self.sel_x, self.sel_y, self.used_y, self.tours
    
    def var_eva_x(self):
        return self.i, self.bt, self.sel_x, self.used_x, self.tours
    
    def var_eva_y(self):
        return self.i, self.bt, self.sel_x, self.sel_y, self.used_y, self.tours


def next_x(iter: Iteration) -> Optional[tuple[Edge, list[int], Edge]]:
    i, bt, sel_y, used_x, tours = iter.var_next_x()
    if i == 0:
        t_start = bt.tbase
        cur_tgraph = bt.tgraph
    else:
        yim1 = sel_y[i-1]
        t_start = yim1[1]
        cur_tgraph, _ = tours[i-1]

    base_nghbs = bt.d_tour[t_start]
    for nghb in base_nghbs:
        maybe_x = (t_start, nghb)

        # condition (a)
        if i <= 1 and edge_in_edge_list(maybe_x, used_x.get(i, [])):
            continue

        if i >= 1:
            # condition (b)
            joined = (bt.tbase, nghb)
            is_tour, maybe_tour = try_is_tour(bt.n, cur_tgraph, maybe_x, joined)
            if not is_tour:
                continue

            # condition (c)
            if edge_in_edge_list(maybe_x, sel_y.values()):
                continue
        else:
            # dummies
            maybe_tour = bt.tour
            joined = (0, 0)

        return maybe_x, maybe_tour, joined

    return None
        
    

def next_y(iter: Iteration) -> Optional[tuple[Edge, float, XResult]]:
    i, bt, sel_x, sel_y, used_y, tours = iter.var_next_y()
    if i == 0:
        partial_gain_im1 = 0
        cur_tgraph = bt.tgraph
    else:
        cur_tgraph, partial_gain_im1 = tours[i-1]
    
    xi = sel_x[i]
    y_start = xi[1]
    y_start_nghbs = bt.d_tour[y_start]

    for node in range(bt.n):
        maybe_yi = (y_start, node)

        # needs to be an edge outside of tour
        if node == y_start or node in y_start_nghbs:
            continue

        # condition (a)
        if i <= 1 and edge_in_edge_list(maybe_yi, used_y.get(i, [])):
            continue
        
        # condition (b)
        gain_i = exchange_gain(xi, maybe_yi, bt.costs)
        if partial_gain_im1 + gain_i <= 0:
            continue

        if i >= 1:
            # condition (c)
            if edge_in_edge_list(maybe_yi, sel_x.values()):
                continue

            # condition (d)
            sel_y[i] = maybe_yi
            yi_tgraph = exchange_edges(cur_tgraph, xi, maybe_yi)
            tours[i] = (yi_tgraph, -1)
            maybe_xi1_res = next_x(iter.with_i(i+1))
            if maybe_xi1_res is None:
                continue
        else:
            maybe_xi1_res = None

        return maybe_yi, gain_i, maybe_xi1_res

    return None


def add_create_edge(i: int, dct: dict[int, set[Edge]], e: Edge):
    if i in dct:
        dct[i].add(e)
    else:
        dct[i] = {e}


# Apologies for the below code, originally it was written with recursion because on small instances it seemed to barely need it
# However, for larger instances it still suffers from Python's recursion limit/lack of recursion optimization, so it was rewritten with as little change as possible

# XEvaluation = tuple[Literal[True], tuple[int, UsedEdges, UsedEdges, SelectedEdges, SelectedEdges, TourState, Optional[XResult]]]
# YEvaluation = tuple[Literal[False], tuple[int, UsedEdges, UsedEdges, SelectedEdges, SelectedEdges, TourState]]
EvalImprovement = tuple[Literal[None], tuple[list[int], float]]


def evaluate_x(iter: Iteration) -> Optional[Iteration]:
    i, bt, sel_x, used_x, tours = iter.var_eva_x()
    xi_result = next_x(iter) if iter.already_xi is None else iter.already_xi
    if xi_result is not None:
        xi, new_tour, joined_i = xi_result

        # check bettter tour!
        if i >= 1:
            if i == 1:
                partial_gain_im1 = 0
            else:
                _, partial_gain_im1 = tours[i-1]

            join_gain_i = exchange_gain(xi, joined_i, bt.costs)

            if partial_gain_im1 + join_gain_i > 0:
                iter.result = (new_tour, join_gain_i)
                return iter

        sel_x[i] = xi
        add_create_edge(i, used_x, xi)
        return iter.with_nxt('y')
    
    if i == 1:
        used_x.pop(1, None)
        # used_y.pop(0, None)
        return iter.with_i(0).with_nxt('y')
    elif i == 0:
        # exhausted all x_0 so go to different t_base
        return None
    else:
        raise ValueError("x_i should exist for higher i!")
    

def evaluate_y(iter: Iteration) -> Iteration:
    i, bt, sel_x, sel_y, used_y, tours = iter.var_eva_y()
    yi_result = next_y(iter)
    # print(f"after next y nxt={iter.next_sel} i={iter.i} used_x={iter.used_x} used_y={iter.used_y} sel_x={iter.sel_x} sel_y={iter.sel_y}")
    
    if yi_result is not None:
        xi = sel_x[i]
        yi, gain_i, xi1_res = yi_result

        sel_y[i] = yi
        add_create_edge(i, used_y, yi)
        
        if i == 0:
            tgraph_im1 = bt.tgraph
            partial_gain_im1 = 0
        else:
            tgraph_im1, partial_gain_im1 = tours[i - 1]

        tgraph_i = exchange_edges(tgraph_im1, xi, yi)
        partial_gain_i = partial_gain_im1 + gain_i
        tours[i] = (tgraph_i, partial_gain_i)

        iter.already_xi = xi1_res
        # print(f"before x i+1 nxt={iter.next_sel} i={iter.i} used_x={iter.used_x} used_y={iter.used_y} sel_x={iter.sel_x} sel_y={iter.sel_y}")
        new_iter = iter.with_i(i+1).with_nxt('x')
        # print(f"incr i with nxt=x i={iter.i} used_x={iter.used_x} used_y={iter.used_y} sel_x={iter.sel_x} sel_y={iter.sel_y}")
        return new_iter
    
    if i >= 2:
        # since we just tried for y_1 we know no more exist
        return iter.with_i(1).with_nxt('x')
    elif i == 1:
        used_y.pop(1, None)
        return iter.with_i(1).with_nxt('x')
    else:
        used_y.pop(0, None)
        return iter.with_i(0).with_nxt('x')
    

def improve_tour(n: int, tbase: int, tour: list[int], costs: list[list[float]]) -> Optional[TourResult]:
    d_tour = dict_tour(tour)
    tgraph = d_tour_to_node_adj_dict(n, d_tour)
    bt = BaseTour(n, tbase, tour, d_tour, tgraph, costs)

    # Very simple iterative implementation of the recursive algorithm, apologies for the ugly code
    base_iter = Iteration(0, 'x', None, dict(), dict(), dict(), dict(), dict(), bt, None)
    queue: list[Iteration] = [base_iter]

    k = 0
    while len(queue) > 0:
        print(f"k {k}")
        k += 1
        current_iter = queue.pop()
        print(f"cur nxt={current_iter.next_sel} i={current_iter.i} used_x={current_iter.used_x} used_y={current_iter.used_y} sel_x={current_iter.sel_x} sel_y={current_iter.sel_y} tours={current_iter.tours}")
        # print(f"i {current_iter.i}")
        if current_iter.next_sel == 'x':
            eval_res = evaluate_x(current_iter)
            # print(eval_res)

            if eval_res is None:
                return None
            elif eval_res.result is not None:
                # print(f"new tour: {eval_res.result[0]}")
                return eval_res.result
            else:
                queue.append(eval_res)

        else:
            eval_res = evaluate_y(current_iter)
            queue.append(eval_res)

    return 



def lin_kernighan(costs: list[list[float]]):
    n = len(costs)
    tour = random_tour(n)
    tour = list(range(n))

    last_tour = tour
    new_tour = tour
    tour_gain = 0

    while new_tour is not None:
        untried_tbase = shuffle_iter(new_tour)
        untried_tbase = new_tour.copy()
        last_tour = new_tour
        new_tour = None

        while len(untried_tbase) > 0:
            tbase = untried_tbase.pop()
            improve_result = improve_tour(n, tbase, last_tour, costs)

            if improve_result is None:
                continue
            else:
                improved_tour, improvement = improve_result
                new_tour = improved_tour
                tour_gain += improvement
                break

    return last_tour, tour_gain


def run():
    # inst_path = get_inst_path()
    inst_path = Path('tsp/burma14.dat')
    graph_l = parse_as_adj_matrix(inst_path)

    print("Value of heuristic")

    # # # Greedy Heuristic: Nearest neighbor
    start_time = perf_counter_ns()
    greedy_tour = nearest_neighbor(graph_l)
    greedy_time = (perf_counter_ns() - start_time) / 1e9
    greedy_length = length_tour(graph_l, greedy_tour)
    print(f"- Nearest Neighbor: {greedy_length}, {greedy_time}s")

    # Lin-Kernighan
    graph_l_fl = parse_as_adj_matrix(inst_path, as_float=True)
    start_time = perf_counter_ns()
    lin_tour, gain = lin_kernighan(graph_l_fl)
    lin_time = (perf_counter_ns() - start_time) / 1e9
    lin_length = length_tour(graph_l, lin_tour)
    print(f"- Lin-Kernighan: {lin_length}, {lin_time}s")
    
    # # # Simulated Annealing1
    # start_time = perf_counter_ns()
    # annealing_tour = simulated_annealing(graph_l, T=100, r=0.95, L=1000)
    # annealing_time = (perf_counter_ns() - start_time) / 1e9
    # annealing_length = length_tour(graph_l, annealing_tour)
    # print(f"- Simulated Annealing: {annealing_length}, {annealing_time}s  (L=1000)")
    
    # # Simulated Annealing2
    # start_time = perf_counter_ns()
    # annealing_tour = simulated_annealing(
    #     graph_l, T=100, r=0.95, L=10000, max_no_improvement=100
    # )
    # annealing_time = (perf_counter_ns() - start_time) / 1e9
    # annealing_length = length_tour(graph_l, annealing_tour)
    # print(
    #     f"- Simulated Annealing: {annealing_length}, {annealing_time}s  (L=10000, max_no_improvement=100)"
    # )

if __name__ == "__main__":
    run()