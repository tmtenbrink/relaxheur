from typing import Iterable, Optional, Union, Literal, Optional
from dataclasses import dataclass
from enum import Enum, auto
import math
from pathlib import Path
import argparse
from heapq import heappop, heappush, heapify
from functools import total_ordering
from random import randrange, random
from time import perf_counter_ns

import gurobipy as gp

from lin_kernighan_2 import lin_kernighan2, normalize_tour, shuffle_normalized


# ==================================================================
# ============================ Parsing =============================
# ==================================================================


Edge = tuple[int, int]
Costs = list[list[float]]

EdgeValue = tuple[float, list[tuple[int, int]]]

# ==================================================================
# ======= Linear Programming Relaxation (Min-Cut Separation) =======
# ==================================================================


def get_edge_names(n: int):
    names = {}
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            names[index] = f"({i}, {j})"

            index += 1

    return names


def get_edges_by_index(n: int) -> dict[int, tuple[int, int]]:
    edges_by_index: dict[int, tuple[int, int]] = {}
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            edges_by_index[index] = (i, j)
            index += 1

    return edges_by_index


class Formulation(Enum):
    EXTENDED = auto()
    CUTTING_PLANE = auto()


@dataclass
class Problem:
    n: int
    # number of edges
    m_edges: int
    # all edges containing a certain vertex (edges according to our index convention)
    all_cuts: list[list[int]]
    # list of coefficients for char_vec (i.e. [1]*m_edges)
    coeff: list[int]
    edge_costs: list[float]

    def all(self) -> tuple[int, int, list[list[int]], list[int]]:
        return self.n, self.m_edges, self.all_cuts, self.coeff


class UnknownVariableError(Exception):
    pass


class TSPModel:
    formulation: Formulation
    model: gp.Model
    p: Problem

    def __init__(
        self, model: gp.Model, n: int, edge_costs: EdgeCosts, formulation: Formulation
    ):
        self.formulation = formulation
        self.model = model
        m_edges = m_edges = (n * (n - 1)) // 2
        all_cuts = edge_idxs_for_all_v(n)
        self.p = Problem(n, m_edges, all_cuts, [1] * len(all_cuts[0]), edge_costs)

    def optimize_with_val(self):
        """RAISES InfeasibleRelaxation when unable to optimize."""
        if self.formulation == Formulation.CUTTING_PLANE:
            self.model = optimize_cut_model(self)
        else:
            self.model.optimize()
        return self.model.ObjVal

    def copy(self):
        return TSPModel(
            self.model.copy(), self.p.n, self.p.edge_costs, self.formulation
        )

    def copy_fix(self, edge_idx: int, fix_val: Literal[0, 1]):
        new_model = self.model.copy()
        char_e = new_model.getVarByName(f"char_{edge_idx}")
        new_model.addConstr(char_e == fix_val, name=f"x_{edge_idx}=={fix_val}")
        return TSPModel(new_model, self.p.n, self.p.edge_costs, self.formulation)

    def char_vec(self) -> list[float]:
        return self.char_vec_values()[1]

    def char_vec_values(self) -> tuple[list[gp.Var], list[float]]:
        m_edges = self.p.m_edges
        var_list = []
        x_values = [0.0] * m_edges
        for e in range(m_edges):
            char_e = self.model.getVarByName(f"char_{e}")
            if char_e is None:
                raise UnknownVariableError(f"Could not find variable char_{e}!")
            var_list.append(char_e)
            x_values[e] = char_e.X
        return var_list, x_values

    def get_tour(self, edge_by_index: Optional[dict[int, tuple[int, int]]] = None):
        """This only works when it is certain it is actually a tour!"""
        _, x_values = self.char_vec_values()
        if edge_by_index is None:
            edge_by_index = get_edges_by_index(self.p.n)

        epsilon = 0.0000001
        tour = []
        queue = []
        for e_i in range(len(x_values)):
            x_val = x_values[e_i]
            if abs(x_val) > epsilon:
                node_1, node_2 = edge_by_index[e_i]
                queue.append((node_1, node_2))

        while queue:
            # 2 cases
            # - node_1 is already in our tour
            #   in that case, since each edge is only once in our list, either:
            #   + node_1 is at the start
            #   + node_1 is at the end
            #   (it cannot be in the middle because then it already has 2 edges connected and it's an invalid tour)
            # - node_1 is not in the tour
            #   in that case, we do same test for node_2
            #   + node_2 is in the tour
            #     * node_2 is at the start
            #     * node_2 is at the end
            #   + node_2 is not in the tour
            #     we put this edge at the end of the queue and continue until we found one that is there
            node_1, node_2 = queue.pop(0)

            if len(tour) == 0:
                tour += [node_1, node_2]
            elif tour[0] == node_1:
                tour.insert(0, node_2)
            elif tour[-1] == node_1:
                tour.append(node_2)
            elif tour[0] == node_2:
                tour.insert(0, node_1)
            elif tour[-1] == node_2:
                tour.append(node_1)
            else:
                queue.append([node_1, node_2])

        # the start/endpoint is going to be added twice, so we have to remove it
        return normalize_tour(tour[1:])

    def print_sol(self):
        char_vec, x_values = self.char_vec_values()
        edge_names = get_edge_names(self.p.n)
        edge_costs = self.p.edge_costs

        epsilon = 0.0000001
        print("Path:")
        cost = 0
        for e_i in range(len(x_values)):
            x_val = x_values[e_i]
            if abs(x_val) > epsilon:
                add_value = f" with value {x_val}"
                print(f"\tedge {edge_names[e_i]} in solution{add_value}")
                cost += edge_costs[e_i] * x_val

        print(f"Cost: {cost}")
        return char_vec, x_values


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


class InfeasibleRelaxation(Exception):
    pass


def cutting_plane_model(n: int, edge_costs: EdgeCosts):
    model = gp.Model("cutting plane")
    # set to dual simplex
    model.setParam("Method", 1)
    # don't log
    model.setParam("LogToConsole", 0)
    print("Setting up cutting plane model...")
    start_build = perf_counter_ns()

    m_edges = (n * (n - 1)) // 2
    char_vec = [model.addVar(lb=0, name=f"char_{e}") for e in range(m_edges)]
    z = gp.LinExpr(edge_costs, char_vec)
    model.setObjective(z, gp.GRB.MINIMIZE)

    model.update()
    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return TSPModel(model, n, edge_costs, Formulation.CUTTING_PLANE)


def get_arc_idxs_for_v(vertex: int, n: int, incoming: bool):
    lower_add = 0 if incoming else 1
    higher_add = 1 - lower_add

    lower_edges = [
        2 * edge_idx(other_v, vertex, n) + lower_add for other_v in range(0, vertex)
    ]
    higher_edges = [
        2 * edge_idx(vertex, other_v, n) + higher_add
        for other_v in range(vertex + 1, n)
    ]

    return lower_edges + higher_edges


def get_in_arc_idxs_for_v(vertex: int, n: int):
    """Return the indexes of the incoming arcs for a vertex, i.e. where v is the head (so second vertex in the ordered pair)."""
    return get_arc_idxs_for_v(vertex, n, True)


def get_out_arc_idxs_for_v(vertex: int, n: int):
    """Return the indexes of the incoming arcs for a vertex, i.e. where v is the head (so second vertex in the ordered pair)."""
    return get_arc_idxs_for_v(vertex, n, False)


def extended_formulation(n: int, edge_costs: EdgeCosts, integer: bool = False):
    model = gp.Model("extended")
    # don't log
    model.setParam("LogToConsole", 0)
    # set to dual simplex
    model.setParam("Method", 1)
    print("Building extended formulation...")
    start_build = perf_counter_ns()

    m_edges = (n * (n - 1)) // 2
    # for i<j (i, j) is the (n-1 + n-2 + ... + n - i) + (j - i -1)th element
    # i.e. index is i(2n-i-1)/2 + (j-i-1)
    edge_vtype = gp.GRB.BINARY if integer else gp.GRB.CONTINUOUS
    char_vec = [
        model.addVar(lb=0, ub=1, name=f"char_{e}", vtype=edge_vtype)
        for e in range(m_edges)
    ]
    num_arcs = 2 * m_edges
    # rows are vertices without r
    # columns are the arcs, with same indexing as above
    # but where (i, j) and (j, i) follow each other (so edge_index * 2 for (i, j)
    # with still i<j
    flow: list[list[gp.Var]] = [
        [model.addVar(lb=0, name=f"f_{s}({a})") for a in range(num_arcs)]
        for s in range(1, n)
    ]
    z = gp.LinExpr(edge_costs, char_vec)
    model.setObjective(z, gp.GRB.MINIMIZE)

    # each column is all the vertices other than the vertex corresponding to column
    # so (n-1), n array
    cut_arr = edge_idxs_for_all_v(n)
    coeff_1 = [1] * len(cut_arr[0])
    for v in range(n):
        v_cut = cut_arr[v]
        cut_x = [char_vec[e] for e in v_cut]
        model.addConstr(gp.LinExpr(coeff_1, cut_x) == 2, name=f"x(delta({v}))==2")
    # we take transpose so we can easily access the rows
    cuts_in = [get_in_arc_idxs_for_v(v, n) for v in range(n)]
    cuts_out = [get_out_arc_idxs_for_v(v, n) for v in range(n)]

    cut_in_r = cuts_in[0]
    cut_out_r = cuts_out[0]
    coeff_r_1_in = [1] * len(cut_in_r)
    coeff_r_1_out = [1] * len(cut_out_r)

    for s in range(1, n):
        f_s = flow[s - 1]
        fs_r_out = [f_s[a] for a in cut_out_r]
        fs_r_in = [f_s[a] for a in cut_in_r]
        out_sum = gp.LinExpr(coeff_r_1_out, fs_r_out)
        in_sum = gp.LinExpr(coeff_r_1_in, fs_r_in)
        model.addConstr(
            out_sum - in_sum >= 2,
            name=f"f_{s}(delta_out(r))-f_{s}(delta_in(r))>=2",
        )

    for s in range(1, n):
        f_s = flow[s - 1]
        # we make sure we get arrays correspond to the right index of the char. vector
        for e in range(m_edges):
            a1 = 2 * e
            a2 = 2 * e + 1
            model.addConstr(f_s[a1] - char_vec[e] <= 0, name=f"f_{a1}<=x({a1})")
            model.addConstr(f_s[a2] - char_vec[e] <= 0, name=f"f_{a2}<=x({a2})")

        for other_v in range(0, n):
            # skip r and s
            if other_v == 0 or other_v == s:
                continue

            cut_out_v = cuts_out[other_v]
            cut_in_v = cuts_in[other_v]
            coeff_v_1_in = [1] * len(cut_in_v)
            coeff_v_1_out = [1] * len(cut_out_v)
            fs_v_out = [f_s[v] for v in cut_out_v]
            fs_v_in = [f_s[v] for v in cut_in_v]

            out_sum = gp.LinExpr(coeff_v_1_out, fs_v_out)
            in_sum = gp.LinExpr(coeff_v_1_in, fs_v_in)

            model.addConstr(
                out_sum - in_sum == 0,
                name=f"f_{s}(delta_out({other_v}))-f_{s}(delta_in({other_v}))==0",
            )
    # We run update to make sure the timer works
    model.update()
    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return TSPModel(model, n, edge_costs, Formulation.EXTENDED)


def build_model(n: int, edge_costs: EdgeCosts, formulation: Formulation, integer=False):
    if formulation == Formulation.EXTENDED:
        return extended_formulation(n, edge_costs, integer)
    elif formulation == Formulation.CUTTING_PLANE:
        if integer:
            raise ValueError("Cutting plane model cannot compute integer value!")
        return cutting_plane_model(n, edge_costs)


# ==================================================================
# =========================== Heuristics ===========================
# ==================================================================


def dict_tour(tour: list[int]) -> dict[int, tuple[int, int]]:
    # adjacent edges in the tour
    d: dict[int, tuple[int, int]] = dict()
    final_i = len(tour) - 1
    for i in range(len(tour)):
        if i == 0:
            other_i1 = 1
            other_i2 = final_i
        elif i == final_i:
            other_i1 = 0
            other_i2 = final_i - 1
        else:
            other_i1 = i + 1
            other_i2 = i - 1

        d[tour[i]] = (tour[other_i1], tour[other_i2])

    return d


def random_tour(n: int) -> list[int]:
    return shuffle_iter(range(n))


def shuffle_iter(itrble: Iterable[int]) -> list[int]:
    """Returns a shuffled copy of the iterable as a list."""
    randoms = [(random(), i) for i in itrble]
    randoms.sort(key=lambda l: l[0])
    shuffled = list(map(lambda r: r[1], randoms))
    return shuffled


Edge = tuple[int, int]
VertTourNghbs = dict[int, tuple[int, int]]
VertNghbs = dict[int, set[int]]

UsedEdges = dict[int, set[Edge]]
SelectedEdges = dict[int, Edge]

# current_graph, partial gain, current tour
IterState = tuple[VertNghbs, float, list[tuple[Edge, Edge]]]


def pick_x0(d_tour: VertTourNghbs, t0: int, costs: list[list[float]]):
    t0_nghbs = d_tour[t0]
    # first one is always t0
    return (t0, t0_nghbs[0]), (t0, t0_nghbs[1])


def exchange_gain(x_orig: Edge, y_repl: Edge, costs: list[list[float]]):
    orig_cost = costs[x_orig[0]][x_orig[1]]
    new_cost = costs[y_repl[0]][y_repl[1]]

    return orig_cost - new_cost


def pick_y0(
    n: int, d_tour: VertTourNghbs, t0: int, t1: int, costs: list[list[float]]
) -> tuple[Optional[Edge], list[Edge], float]:
    in_tour = d_tour[t1]
    other_tour_nghb = in_tour[0] if in_tour[0] != t0 else in_tour[1]
    possible_y0 = [
        (t1, i) for i in range(n) if i != t0 and i != t1 and i != other_tour_nghb
    ]

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


def edge_eq(e1: Edge, e2: Edge) -> bool:
    return e1 == e2 or (e1[1], e1[0]) == e2


def edge_in_edge_list(e: Edge, edge_list: Iterable[Edge]) -> bool:
    return e in edge_list or (e[1], e[0]) in edge_list


@dataclass
class TourState:
    n: int
    tbase: int
    tour: list[int]
    d_tour: VertTourNghbs
    tgraph: VertNghbs
    costs: list[list[float]]
    iter: dict[int, IterState]


XResult = Optional[tuple[Edge, list[int], Edge]]


def next_x(
    i: int,
    used_x: UsedEdges,
    sel_y: SelectedEdges,
    ts: TourState,
    tgraph_moves: Optional[tuple[VertNghbs, list[tuple[Edge, Edge]]]] = None,
) -> Optional[tuple[Edge, list[int], Edge]]:
    if i == 0:
        t_start = ts.tbase
        cur_tgraph = ts.tgraph
        moves = []
    else:
        yim1 = sel_y[i - 1]
        t_start = yim1[1]
        if tgraph_moves is None:
            cur_tgraph, _, moves = ts.iter[i - 1]
        else:
            cur_tgraph, moves = tgraph_moves

    base_nghbs = ts.d_tour[t_start]
    for nghb in base_nghbs:
        maybe_x = (t_start, nghb)

        # condition (a)
        if i <= 1 and edge_in_edge_list(maybe_x, used_x.get(i, [])):
            continue

        if i >= 1:
            # condition (b)
            joined = (ts.tbase, nghb)

            # this is not possible of course
            if joined[0] == joined[1]:
                continue

            # replacing it with the same edge doesn't work of course
            if joined == maybe_x:
                continue

            is_tour, maybe_tour = try_is_tour(
                ts.n, cur_tgraph, maybe_x, joined, moves, ts.tgraph
            )
            if not is_tour:
                continue

            # condition (c)
            if edge_in_edge_list(maybe_x, sel_y.values()):
                continue
        else:
            # dummies
            maybe_tour = ts.tour
            joined = (0, 0)

        return maybe_x, maybe_tour, joined

    return None


def next_y(
    i: int,
    used_x: UsedEdges,
    used_y: UsedEdges,
    sel_x: SelectedEdges,
    sel_y: SelectedEdges,
    ts: TourState,
) -> Optional[tuple[Edge, float, XResult]]:
    if i == 0:
        partial_gain_im1 = 0
        cur_tgraph = ts.tgraph
        moves = []
    else:
        cur_tgraph, partial_gain_im1, moves = ts.iter[i - 1]

    xi = sel_x[i]
    y_start = xi[1]
    y_start_nghbs = ts.d_tour[y_start]

    for node in range(ts.n):
        maybe_yi = (y_start, node)

        # needs to be an edge outside of tour
        if node == y_start or node in y_start_nghbs:
            continue

        # condition (a)
        if i <= 1 and edge_in_edge_list(maybe_yi, used_y.get(i, [])):
            continue

        # condition (b)
        gain_i = exchange_gain(xi, maybe_yi, ts.costs)
        if partial_gain_im1 + gain_i <= 0:
            continue

        if i >= 1:
            # condition (c)
            if edge_in_edge_list(maybe_yi, sel_x.values()):
                continue

            # condition (d)
            sel_y[i] = maybe_yi
            yi_tgraph = exchange_edges(cur_tgraph, xi, maybe_yi)

            maybe_xi1_res = next_x(
                i + 1,
                used_x,
                sel_y,
                ts,
                tgraph_moves=(yi_tgraph, moves + [(xi, maybe_yi)]),
            )
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


XEvaluation = tuple[
    Literal[True],
    tuple[
        int,
        UsedEdges,
        UsedEdges,
        SelectedEdges,
        SelectedEdges,
        TourState,
        Optional[XResult],
    ],
]
YEvaluation = tuple[
    Literal[False],
    tuple[int, UsedEdges, UsedEdges, SelectedEdges, SelectedEdges, TourState],
]
EvalImprovement = tuple[Literal[None], tuple[list[int], float]]


def try_is_tour(
    n: int,
    tour_n_a_dict: VertNghbs,
    x_remove: Edge,
    y_repl: Edge,
    moves: list[tuple[Edge, Edge]],
    base_tgraph: VertNghbs,
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


def is_node_edges_tour(
    n: int, node_edges: VertNghbs, x_remove: Edge, y_repl: Edge
) -> tuple[bool, list[int]]:
    maybe_tour_n_a_dict = exchange_edges(node_edges, x_remove, y_repl, copy=False)
    tour = vert_ngbhs_to_tour_seen(maybe_tour_n_a_dict)
    exchange_edges(maybe_tour_n_a_dict, y_repl, x_remove, copy=False)

    if len(tour) != n:
        return False, []
    return True, tour


def copy_node_edges(node_edges: VertNghbs) -> VertNghbs:
    node_edge_copy = dict()
    for t in node_edges:
        node_edge_copy[t] = node_edges[t].copy()

    return node_edge_copy


def evaluate_x(
    i: int,
    used_x: UsedEdges,
    used_y: UsedEdges,
    sel_x: SelectedEdges,
    sel_y: SelectedEdges,
    ts: TourState,
    xi_from_y: Optional[XResult] = None,
) -> Optional[Union[YEvaluation, EvalImprovement]]:
    xi_result = next_x(i, used_x, sel_y, ts) if xi_from_y is None else xi_from_y
    if xi_result is not None:
        xi, new_tour, joined_i = xi_result

        # check bettter tour!
        if i >= 1:
            if i == 1:
                partial_gain_im1 = 0
            else:
                _, partial_gain_im1, _ = ts.iter[i - 1]

            join_gain_i = exchange_gain(xi, joined_i, ts.costs)

            if partial_gain_im1 + join_gain_i > 0:
                return (None, (new_tour, join_gain_i))

        sel_x[i] = xi
        add_create_edge(i, used_x, xi)
        return (False, (i, used_x, used_y, sel_x, sel_y, ts))

    if i == 1:
        used_x.pop(1, None)
        # used_y.pop(0, None)
        return (False, (0, used_x, used_y, sel_x, sel_y, ts))
    elif i == 0:
        # exhausted all x_0 so go to different t_base
        return None
    else:
        raise ValueError("x_i should exist for higher i!")


def evaluate_y(
    i: int,
    used_x: UsedEdges,
    used_y: UsedEdges,
    sel_x: SelectedEdges,
    sel_y: SelectedEdges,
    ts: TourState,
) -> XEvaluation:
    yi_result = next_y(i, used_x, used_y, sel_x, sel_y, ts)

    if yi_result is not None:
        xi = sel_x[i]
        yi, gain_i, xi1_res = yi_result

        sel_y[i] = yi
        add_create_edge(i, used_y, yi)

        if i == 0:
            tgraph_im1 = ts.tgraph
            partial_gain_im1 = 0
            moves = []
        else:
            tgraph_im1, partial_gain_im1, moves = ts.iter[i - 1]

        tgraph_i = exchange_edges(tgraph_im1, xi, yi)
        partial_gain_i = partial_gain_im1 + gain_i
        ts.iter[i] = (tgraph_i, partial_gain_i, moves + [(xi, yi)])

        i = i + 1

        return (True, (i, used_x, used_y, sel_x, sel_y, ts, xi1_res))

    if i >= 2:
        # since we just tried for y_1 we know no more exist
        return (True, (1, used_x, used_y, sel_x, sel_y, ts, None))
    elif i == 1:
        used_y.pop(1, None)
        return (True, (1, used_x, used_y, sel_x, sel_y, ts, None))
    else:
        used_y.pop(0, None)
        return (True, (0, used_x, used_y, sel_x, sel_y, ts, None))


def improve_tour(n: int, tbase: int, tour: list[int], costs: list[list[float]]):
    d_tour = dict_tour(tour)
    tgraph = d_tour_to_node_adj_dict(n, d_tour)
    ts = TourState(n, tbase, tour, d_tour, tgraph, costs, dict())

    # Very simple iterative implementation of the recursive algorithm, apologies for the ugly code
    queue: list[Union[XEvaluation, YEvaluation]] = [
        (True, (0, dict(), dict(), dict(), dict(), ts, None))
    ]

    while len(queue) > 0:
        eval_q = queue.pop()

        if eval_q[0] is True:
            _, x_eval = eval_q
            eval_res = evaluate_x(*x_eval)

            if eval_res is None:
                return None
            elif eval_res[0] is None:
                _, improvement = eval_res
                return improvement
            else:
                queue.append(eval_res)

        else:
            _, y_eval = eval_q
            eval_res = evaluate_y(*y_eval)
            queue.append(eval_res)

    return


def normalize_tour_repr(tour: list[int]):
    max_size = len(tour)

    bytes_per_n = (max_size.bit_length() + 7) // 8

    i_zero = tour.index(0)
    prev_part = tour[:i_zero]
    normalized_tour = tour[i_zero:] + prev_part
    byte_seq = [i.to_bytes(bytes_per_n, byteorder="big") for i in normalized_tour]
    tour_bytes = b"".join(byte_seq)
    tour_bytes_rev = byte_seq[0] + b"".join(reversed(byte_seq[1:]))

    return min(tour_bytes, tour_bytes_rev)


def lin_kernighan(
    costs: list[list[float]], start_tour: Optional[list[int]] = None, no_random=False
):
    n = len(costs)

    if start_tour is not None:
        tour = start_tour
    elif no_random:
        tour = list(range(n))
    else:
        tour = random_tour(n)

    last_tour = tour
    new_tour = tour
    tour_gain = 0

    while new_tour is not None:
        if no_random:
            untried_tbase = new_tour.copy()
        else:
            untried_tbase = shuffle_iter(new_tour)
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


def check_fixed_tour(tour: list[int], fixed_one: list[Edge], fixed_zero: list[Edge]):
    tour_edges_l = []
    for i in range(len(tour)):
        if i == 0:
            continue
        tour_edges_l.append((tour[i], tour[i - 1]))
        tour_edges_l.append((tour[i - 1], tour[i]))
    tour_edges_l.append((tour[0], tour[-1]))
    tour_edges_l.append((tour[-1], tour[0]))
    tour_edges = set(tour_edges_l)

    for e in fixed_zero:
        if e in tour_edges:
            raise ValueError(f"Edge {e} is not allowed in tour but was found!")
    for e in fixed_one:
        if e not in tour_edges:
            raise ValueError(f"Edge {e} is required in tour but was not found!")


# ==================================================================
# ======================== Branch and Bound ========================
# ==================================================================


@total_ordering
@dataclass
class Subproblem:
    fixed_one: list[Edge]
    fixed_zero: list[Edge]
    parent_lb: float
    branch_e_idx: int
    # upward if the branch var has been set to one, downward otherwise
    upward_branch: bool
    parent_branch_e_val: float
    model: TSPModel
    # with fixed edges set to inf/ninf
    heur_costs: list[list[float]]

    # we define lt method so that they can be compared by heapsort
    def __lt__(self, other):
        return self.parent_lb < other.parent_lb


# LP time, heur time


def gurobi_integer(inst: list[list[float]]):
    """Runs Gurobi to solve TSP
    ins: adjacency matrix"""

    n = len(inst)
    edge_costs = compute_edge_costs(inst)
    ext_model = build_model(n, edge_costs, Formulation.EXTENDED, integer=True)
    print(f"Optimizing integer extended formulation...")
    start_opt = perf_counter_ns()
    value = ext_model.optimize_with_val()
    opt_time = (perf_counter_ns() - start_opt) / 1e9
    tour = ext_model.get_tour()
    print(f"\t- integer optimal: {value}")
    print(f"\t- optimal tour: {tour}")
    print(f"Optimizing extended integer model using only Gurobi took {opt_time} s.")


def feasible_tour(
    costs: list[list[float]],
    fixed_one: list[tuple[int, int]],
    fixed_zero: list[tuple[int, int]],
    starting_node=None,
):
    if starting_node is None:
        starting_node = randrange(0, len(costs), 1)

    # Count the number of times a node is in fixed_one
    count_dict = {}
    for tup in fixed_one:
        for num in tup:
            count_dict[num] = count_dict.get(num, 0) + 1

    if any(count > 2 for count in count_dict.values()):
        raise ValueError(
            "Error: problem infeasible. Some nodes have more than 2 fixed edges."
        )

    # Nodes that appear twice in fixed_ones can only be added using fixed edges, not using nearest neighbor
    nodes_appear_twice = [node for node, count in count_dict.items() if count == 2]

    n = len(costs)
    remaining_nodes = list(range(n))

    cur_node = starting_node
    remaining_nodes.remove(cur_node)
    tour = [cur_node]

    for _ in range(n - 1):
        # if fixed -> add that one
        cur_node = tour[-1]

        if any(cur_node in edge for edge in fixed_one):
            for edge in fixed_one:
                if cur_node in edge:
                    next_node = edge[0] if cur_node == edge[1] else edge[1]
                    tour.append(next_node)
                    fixed_one.remove(edge)
                    remaining_nodes.remove(next_node)

        # else -> add neareast neighbor
        else:
            # get the costs for all nodes
            costs_row = costs[cur_node]

            # Check which nodes can be visited
            neighbors = [
                node
                for node in remaining_nodes
                if (cur_node, node) not in fixed_zero
                and (node, cur_node) not in fixed_zero
                and node not in nodes_appear_twice
            ]
            if not neighbors:
                break  # not possible to find new node

            # get a start index
            next_node = neighbors[0]

            for n in neighbors[1:]:
                if costs_row[n] < costs_row[next_node]:
                    next_node = n

            tour.append(next_node)
            remaining_nodes.remove(next_node)

    return tour


def main():
    # inst_path = get_inst_path()
    inst_path = Path("tsp/gr96.dat")
    graph_l = parse_as_adj_matrix(inst_path)

    do_branch_and_bound(graph_l)

    gurobi_integer(graph_l)


main()
