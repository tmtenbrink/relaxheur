from dataclasses import dataclass
from pathlib import Path
import argparse
from heapq import heappop, heappush, heapify
from functools import total_ordering
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


def parse_line(ln: str) -> list[float]:
    return list(map(lambda i: float(i), ln.rstrip().split(" ")))


def get_inst_path():
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument("inst_name", help="Filename of instance")
    args = parser.parse_args()

    return Path(args.inst_name)

def parse_as_adj_matrix(inst_path: Path) -> list[list[float]]:
    with open(inst_path, "r") as f:
        lines = f.readlines()

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

Edge = tuple[int, int]

@total_ordering
@dataclass
class Subproblem:
    fixed_one: set[Edge]
    fixed_zero: set[Edge]
    lb: float       

    # we define lt method so that they can be compared by heapsort
    def __lt__(self, other):
        return self.lb < other.lb


class InfeasibleSubproblem(Exception):
    pass


def compute_bounds(inst: list[list[float]], problem: Subproblem):
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


def do_branch_and_bound(inst: list[list[float]]):
    # Global upper and lower bounds, and best solution found so far.
    # inst is adjacency matrix
    sum_ = 0
    # cost is definitely lower than the sum of all the costs
    for lst in inst:
        sum_ += sum(lst)
    global_ub = sum_
    global_lb = -1.0
    best_solution = []

    n = len(inst)

    # ====== AANPASSEN  LB <-> UB

    # Initialization.
    active_nodes: list[Subproblem] = [Subproblem(set(), set(), global_lb)]
    # just to indicate we are working with a list with the heap property
    heapify(active_nodes)

    # Main loop.
    while active_nodes:
        # Select an active node to process. We choose the one with the lowest lower bound (hopefully this will also 
        # have a good upper bound, alllowing us to prune faster)
        problem = heappop(active_nodes)

        # Process the node.
        try:
            lb, ub, tour = compute_bounds(inst, problem)
        except InfeasibleSubproblem:
            # Pruned by infeasibility.
            continue

        # Update global upper bound.
        if ub < global_ub:
            global_ub = ub
            best_solution = tour 
            print("Improved lower bound:", global_ub)

        # by heap property we have that active_nodes[0] has the lowest lower bound
        new_global_lb = min(lb, active_nodes[0].lb)

        if new_global_lb > global_lb:
            global_lb = new_global_lb
            print("Improved lower bound:", global_lb)

        # Prune by bound
        if lb > global_ub:
            continue

        # Prune by optimality
        if lb == ub:
            continue

        # Select edge for split and perform the split.
        # use LP relaxation -> select variable to split, one that is not integer
        # TODO Fix this
        for i in range(n):
            if i not in problem.fixed_one and i not in problem.fixed_zero:
                Pl = (list(p_fixed_one) + [i], p_fixed_zero, lb)
                Pr = (p_fixed_one, list(p_fixed_zero) + [i], lb)
                active_nodes += [Pl, Pr]
        else:
            raise RuntimeError("no variable to fix; this is a bug")

        
    assert global_ub >= global_lb

    # Check that the solution is truly feasible.
    # if infeasible:
    #    raise RuntimeError('solution is infeasible; this is a bug')

    # Return optimal solution.
    return global_ub, best_solution


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
