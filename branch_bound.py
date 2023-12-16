from dataclasses import dataclass
from pathlib import Path
import argparse
from heapq import heappop, heappush, heapify
from functools import total_ordering
from random import randrange

import gurobipy as gp

from typing import Iterable, Optional, Union, Literal, Optional, Union, overload

from heuristics import check_fixed_tour, lin_kernighan, lin_kernighan_fixed
from relax import extended_formulation

from relax_bb import (
    EdgeCosts,
    InfeasibleRelaxation,
    cutting_plane_model,
    Costs,
    TSPModel,
    edge_idx,
    get_edges_by_index
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

Edge = tuple[int, int]

# =====================================
# ========== Code Relaxation ==========
# =====================================

def compute_edge_costs(costs: Costs) -> EdgeCosts:
    n = len(costs)
    m_edges = (n * (n - 1)) // 2
    edge_costs = [0.0]*m_edges

    for row_num, row in enumerate(costs):
        target_start = row_num * (2 * n - row_num - 1) // 2
        target_end = target_start + (n - row_num  - 1)
        for i, e in enumerate(range(target_start, target_end)):
            edge_costs[e] = row[row_num+1+i]

    return edge_costs


# ====================================
# ========== Code Heuristic ==========
# ====================================

def length_tour(graph: list[list[float]], tour: list[int]):
    length = graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        length += graph[tour[i]][tour[i + 1]]

    return length


# ======================================
# ========== Branch and Bound ==========
# ======================================



@total_ordering
@dataclass
class Subproblem:
    fixed_one: list[Edge]
    fixed_zero: list[Edge]
    lb: float
    model: TSPModel
    # with fixed edges set to inf/ninf
    heur_costs: list[list[float]]

    # we define lt method so that they can be compared by heapsort
    def __lt__(self, other):
        return self.lb < other.lb


def compute_bounds(costs: list[list[float]], problem: Subproblem, best_solution: list[int]):
    # Solve relaxation to find LB
    # RAISES InfeasibleRelaxation 
    lb = problem.model.optimize_with_val()

    heur_tour, ub = heuristic(problem, costs, best_solution)

    return lb, ub, heur_tour


def intialize(costs: Costs) -> Subproblem:
    n = len(costs)
    edge_costs = compute_edge_costs(costs)
    cut_model = cutting_plane_model(n, edge_costs)

    return Subproblem([], [], -1, cut_model, copy_costs(costs))


def heuristic(problem: Subproblem, costs: Costs, best_solution: list[int]) -> tuple[list[int], float]:
    tour, _ = lin_kernighan(problem.heur_costs, start_tour=best_solution)
    check_fixed_tour(tour, problem.fixed_one, problem.fixed_zero)
    ub = length_tour(costs, tour)
    return tour, ub


def initial_ub(inst: Costs) -> float:
    ub = 0.0
    # cost is definitely lower than the sum of all the costs
    for lst in inst:
        ub += sum(lst)
    return ub

def copy_costs(costs: Costs):
    return [row.copy() for row in costs]

def heur_fix(heur_costs: Costs, edge: Edge, fix_val: Literal[0, 1]):
    new_heur_costs = copy_costs(heur_costs)
    if fix_val == 1:
        new_heur_costs[edge[0]][edge[1]] = -float('inf')
        new_heur_costs[edge[1]][edge[0]] = -float('inf')
    elif fix_val == 0:
        new_heur_costs[edge[0]][edge[1]] = float('inf')
        new_heur_costs[edge[1]][edge[0]] = float('inf')
    return new_heur_costs

def do_branch_and_bound(inst: Costs):
    # Global upper and lower bounds, and best solution found so far.
    # inst is adjacency matrix
    
    global_problem = intialize(inst)
    n = len(inst)

    global_lb = global_problem.lb
    global_ub = initial_ub(inst)
    
    best_solution = list(range(n))

    edges_by_index = get_edges_by_index(n)

    # Initialization.
    active_nodes: list[Subproblem] = [global_problem]
    # just to indicate we are working with a list with the heap property
    heapify(active_nodes)

    # Main loop.
    while active_nodes:
        # Select an active node to process. We choose the one with the lowest lower bound (hopefully this will also 
        # have a good upper bound, alllowing us to prune faster)
        problem = heappop(active_nodes)

        # Process the node.
        try:
            lb, ub, tour = compute_bounds(inst, problem, best_solution)
            print(f"New suproblem: lb={lb}, ub={ub} (glb={global_lb}, gub={global_ub})")
        except InfeasibleRelaxation:
            print("Infeasible suproblem...")
            # Pruned by infeasibility.
            continue

        # Update global upper bound.
        if ub < global_ub:
            global_ub = ub
            best_solution = tour 
            print("Improved upper bound:", global_ub)

        # by heap property we have that active_nodes[0] has the lowest lower bound
        new_global_lb = min(lb, active_nodes[0].lb) if active_nodes else lb

        # if it's greater than global ub we've already found an optimum and we're just making things worse
        if new_global_lb > global_lb and new_global_lb < global_ub:
            global_lb = new_global_lb
            print("Improved lower bound:", global_lb)

        # Prune by bound
        if lb > global_ub:
            print("Pruning suproblem by worse lower bound...")
            continue
            

        # Prune by optimality
        if lb == ub:
            print("Optimal...")
            continue

        # split into two based on first non-integer edge
        # print(f"current: {problem.model.char_vec_values()}")
        for e_idx, e_val in enumerate(problem.model.char_vec()):
            epsilon = 0.0000001
            if abs(round(e_val) - e_val) > epsilon:
                print(f"edge {e_idx} with value {e_val} is not integer. Splitting...")
                edge = edges_by_index[e_idx]
                new_model_l = problem.model.copy_fix(e_idx, 1)
                new_model_r = problem.model.copy_fix(e_idx, 0)

                fixed_one_l = problem.fixed_one + [edge]
                fixed_zero_r = problem.fixed_zero + [edge]

                heur_costs_l = heur_fix(problem.heur_costs, edge, 1)
                heur_costs_r = heur_fix(problem.heur_costs, edge, 0)

                problem_l = Subproblem(fixed_one_l, problem.fixed_zero, lb, new_model_l, heur_costs_l)
                problem_r = Subproblem(problem.fixed_one, fixed_zero_r, lb, new_model_r, heur_costs_r)
                
                # we ensure they get added in the right spot to keep the heap invariant
                heappush(active_nodes, problem_l)
                heappush(active_nodes, problem_r)
                break
                
        else:
            # we've found an integer solution that the heuristic couldn't find
            print(f"Integer LP solution...")
            if lb < global_ub:
                tour = problem.model.get_tour(edges_by_index)
                global_lb = ub
                best_solution = tour 
                print(f"Improved upper bound {global_ub}.")
            
            # problem.model.print_sol()
            # raise RuntimeError("no variable to fix; this is a bug")

        
    assert global_ub >= global_lb

    # Check that the solution is truly feasible.
    # if infeasible:
    #    raise RuntimeError('solution is infeasible; this is a bug')

    # Return optimal solution.
    return global_ub, best_solution


def gurobi_integer(inst: list[list[float]]):
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
    print("\nOptimal solution =", opt, " tour =", sol)

    print("\nNow we run Gurobi...\n")
    gurobi_integer(graph_l)


main()
