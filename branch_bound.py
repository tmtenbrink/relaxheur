from pathlib import Path
import argparse

from random import randrange

import gurobipy as gp

from typing import Iterable, Optional, Union, Literal, Optional, Union, overload

from heuristics import lin_kernighan

from relax import TSPModel
from relax import (
    compute_edge_costs,
    cutting_plane_model,
    extended_formulation,
)


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

    if as_float:
        adj_matrix = list(map(lambda ln: parse_line(ln, True), lines[1:]))
    else:
        adj_matrix = list(map(lambda ln: parse_line(ln), lines[1:]))

    return adj_matrix


# =====================================
# ========== Code Relaxation ==========
# =====================================


# ====================================
# ========== Code Heuristic ==========
# ====================================
def length_tour(graph, tour):
    length = graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        length += graph[tour[i]][tour[i + 1]]

    return length


# ======================================
# ========== Branch and Bound ==========
# ======================================


class Subproblem:
    def __init__(self, fixed_one, fixed_zero, lb):
        self.fixed_one = set(fixed_one)
        self.fixed_zero = set(fixed_zero)
        self.lb = lb


class InfeasibleSubproblem(Exception):
    pass


def compute_bounds(inst, P):
    # if infeasible:
    #    raise InfeasibleSubproblem()

    # Solve relaxation to find LB
    LP_inst = inst  # copy + fix edges -> constraints

    n = len(inst)
    edge_costs = compute_edge_costs(LP_inst)
    cut_model = cutting_plane_model(n, edge_costs)
    # code aanpassen naar dat hij niet de obj value print, maar alleen obj value geeft
    cut_model.optimize()
    lb = cut_model.model.ObjVal  # werkt dit??

    # Solve heuristic to find UB
    # fix edges -> set high/low value
    heur_inst = inst  # + constraints

    lin_tour, gain = lin_kernighan(heur_inst)
    ub = length_tour(heur_inst, lin_tour)

    # ?? output lin_tour as feasible solution ??
    return lb, ub, lin_tour


def do_branch_and_bound(inst):
    # Global upper and lower bounds, and best solution found so far.
    # inst is adjacency matrix
    sum_ = 0
    for lst in inst:
        sum_ += sum(lst)
    global_ub = sum_
    global_lb = -1
    best_solution = []

    # ====== AANPASSEN  LB <-> UB

    # After processing a node, the new global lower bound is the
    # minimum of the lower bounds of active nodes and integer
    # nodes. Since integer nodes get out of the list of active nodes,
    # we keep the minimum lower bound of integer nodes in the
    # following variable.
    integer_node_bound = -1

    # Initialization.
    active_nodes = [Subproblem([], [], global_ub)]  # ==== checken

    # Main loop.
    while active_nodes:
        # Select an active node to process.
        P = active_nodes[0]
        active_nodes = active_nodes[1:]

        # Process the node.
        try:
            lb, ub, items = compute_bounds(inst, P)  # ====== AANPASSEN items
        except InfeasibleSubproblem:
            # Pruned by infeasibility.
            continue

        # Update global upper bound.
        if ub < global_ub:
            global_ub = ub
            best_solution = list(P.fixed_one) + items  # ====== AANPASSEN
            print("Improved lower bound:", global_ub)

        # Update global lower bound.
        if lb == ub and ub > integer_node_bound:  # ???
            integer_node_bound = lb

        if active_nodes:
            new_global_lb = min(lb, integer_node_bound, min(P.lb for P in active_nodes))
        else:
            new_global_lb = min(lb, integer_node_bound)

        if new_global_lb > global_lb:
            global_lb = new_global_lb
            print("Improved lower bound:", global_lb)

        # Prune by bound?
        if lb > global_ub:
            continue

        # Prune by optimality?
        if lb == ub:
            continue

        # Select variable for split and perform the split.
        # use LP relaxation -> select variable to split, one that is not integer
        for i in range(inst.n):
            if i not in P.fixed_one and i not in P.fixed_zero:
                break
        else:
            raise RuntimeError("no variable to fix; this is a bug")

        Pl = Subproblem(list(P.fixed_one) + [i], P.fixed_zero, lb)
        Pr = Subproblem(P.fixed_one, list(P.fixed_zero) + [i], lb)
        active_nodes += [Pl, Pr]

    assert global_ub >= global_lb

    # Check that the solution is truly feasible.
    # if infeasible:
    #    raise RuntimeError('solution is infeasible; this is a bug')

    # Return optimal solution.
    return ub, best_solution


def run_gurobi(inst):
    """Runs Gurobi to solve TSP
    ins: adjacency matrix"""

    n = len(inst)
    edge_costs = compute_edge_costs(inst)
    ext_model = extended_formulation(n, edge_costs)
    ext_model.optimize()


def main():
    # inst_path = get_inst_path()
    inst_path = Path("tsp/gr48.dat")
    graph_l = parse_as_adj_matrix(inst_path)

    print("Here is the output of the branch-and-bound method")
    opt, sol = do_branch_and_bound(graph_l)
    print("\nOptimal solution =", opt, " items =", sol)

    print("\nNow we run Gurobi...\n")
    try:
        run_gurobi(graph_l)
    except RuntimeError:
        print("error")


main()
