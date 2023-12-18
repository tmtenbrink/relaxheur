from dataclasses import dataclass
import math
from pathlib import Path
import argparse
from heapq import heappop, heappush, heapify
from functools import total_ordering
from random import randrange, random
from time import perf_counter_ns

import gurobipy as gp

from typing import Iterable, Optional, Union, Literal, Optional, Union, overload

from heuristics import check_fixed_tour, lin_kernighan, lin_kernighan_fixed
from pseudocost_branching import PseudoList, pseudocost, update_eta_sigma
from relax import extended_formulation

from relax_bb import (
    EdgeCosts,
    Formulation,
    InfeasibleRelaxation,
    build_model,
    cutting_plane_model,
    Costs,
    TSPModel,
    edge_idx,
    get_edges_by_index,
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
    edge_costs = [0.0] * m_edges

    for row_num, row in enumerate(costs):
        target_start = row_num * (2 * n - row_num - 1) // 2
        target_end = target_start + (n - row_num - 1)
        for i, e in enumerate(range(target_start, target_end)):
            edge_costs[e] = row[row_num + 1 + i]

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
Timer = tuple[int, int]


def compute_bounds(
    costs: list[list[float]],
    problem: Subproblem,
    best_solution: list[int],
    best_solution_cost: float,
    timer: Timer,
    heur_amount_i: int,
):
    # Solve relaxation to find LB
    # RAISES InfeasibleRelaxation
    before = perf_counter_ns()
    lb = problem.model.optimize_with_val()
    after_lp = perf_counter_ns()

    heur_tour, ub = heuristic(problem, costs, best_solution, best_solution_cost, heur_amount_i)
    after_heur = perf_counter_ns()

    timer = (timer[0] + (after_lp - before), timer[1] + (after_heur - after_lp))

    return lb, ub, heur_tour, timer


def intialize(costs: Costs, formulation=Formulation.CUTTING_PLANE) -> Subproblem:
    n = len(costs)
    edge_costs = compute_edge_costs(costs)
    lp_model = build_model(n, edge_costs, formulation)

    return Subproblem([], [], -1, -1, False, -1, lp_model, copy_costs(costs))


def heuristic(
    problem: Subproblem, costs: Costs, best_solution: list[int], best_sol_cost: float, heur_amount: int = 5
) -> tuple[list[int], float]:
    if heur_amount <= 0:
        return best_solution, best_sol_cost
    
    best_heur_tour, _ = lin_kernighan(problem.heur_costs, start_tour=best_solution, no_random=True)
    best_cost = length_tour(costs, best_heur_tour)
    lengths = []
    for _ in range(heur_amount):
        tour, _ = lin_kernighan(problem.heur_costs, start_tour=best_solution)
        tour_cost = length_tour(costs, tour)
        lengths.append(tour_cost)
        if tour_cost < best_cost:
            best_heur_tour = tour
            best_cost = tour_cost
    # print(lengths)
    check_fixed_tour(best_heur_tour, problem.fixed_one, problem.fixed_zero)
    ub = length_tour(costs, best_heur_tour)
    return best_heur_tour, ub


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
        new_heur_costs[edge[0]][edge[1]] = -float("inf")
        new_heur_costs[edge[1]][edge[0]] = -float("inf")
    elif fix_val == 0:
        new_heur_costs[edge[0]][edge[1]] = float("inf")
        new_heur_costs[edge[1]][edge[0]] = float("inf")
    return new_heur_costs


def branch_variable(problem: Subproblem, edge: Edge, e_idx: int, branch_e_val: float, parent_lb: float):
    new_model_l = problem.model.copy_fix(e_idx, 1)
    new_model_r = problem.model.copy_fix(e_idx, 0)

    fixed_one_l = problem.fixed_one + [edge]
    fixed_zero_r = problem.fixed_zero + [edge]

    heur_costs_l = heur_fix(problem.heur_costs, edge, 1)
    heur_costs_r = heur_fix(problem.heur_costs, edge, 0)

    problem_l = Subproblem(
        fixed_one_l, problem.fixed_zero, parent_lb, e_idx, True, branch_e_val, new_model_l, heur_costs_l
    )
    problem_r = Subproblem(
        problem.fixed_one, fixed_zero_r, parent_lb, e_idx, False, branch_e_val, new_model_r, heur_costs_r
    )

    return problem_l, problem_r


def find_branch_variable(
    problem: Subproblem,
    edges_by_index: dict[int, tuple[int, int]],
    pseudo_plus: PseudoList,
    pseudo_min: PseudoList,
    method: Literal["pseudocost", "first-non-integer"] = "pseudocost",
) -> Optional[tuple[Edge, int, float]]:
    if method == "first-non-integer":
        for e_idx, e_val in enumerate(problem.model.char_vec()):
            epsilon = 0.0000001
            if abs(round(e_val) - e_val) > epsilon:
                # print(f"edge {e_idx} with value {e_val} is not integer. Splitting...")
                return edges_by_index[e_idx], e_idx, e_val
    elif method == "pseudocost":
        e_idx, e_val = pseudocost(problem.model.char_vec(), pseudo_plus, pseudo_min)
        if e_idx == -1:
            return None
        return edges_by_index[e_idx], e_idx, e_val

    return None


def do_branch_and_bound(inst: Costs):
    # Global upper and lower bounds, and best solution found so far.
    # inst is adjacency matrix

    global_problem = intialize(inst)
    n = len(inst)

    global_lb = global_problem.parent_lb
    global_ub = initial_ub(inst)

    best_solution = list(range(n))

    edges_by_index = get_edges_by_index(n)

    m_edges = global_problem.model.p.m_edges

    pseudo_plus = ([0.0]*m_edges, [0.0]*m_edges, list[int]())
    pseudo_min = ([0.0]*m_edges, [0.0]*m_edges, list[int]())

    # Initialization.
    active_nodes: list[Subproblem] = [global_problem]
    # we use a priority queue since we use the heuristic to quickly get good upper bounds
    # this is best-first search
    heapify(active_nodes)

    start_opt = perf_counter_ns()
    timer: Timer = (0, 0)

    heur_amount = 25

    # Main loop.
    while active_nodes:
        # Select an active node to process. We choose the one with the lowest lower bound (hopefully this will also
        # have a good upper bound, alllowing us to prune faster)
        problem = heappop(active_nodes)

        heur_amount_i = round(math.ceil(random() < heur_amount)*max(heur_amount, 0.501))
        heur_amount = max(0.25, heur_amount/1.3)

        # Process the node.
        try:
            lb, ub, tour, timer = compute_bounds(inst, problem, best_solution, global_ub, timer, heur_amount_i)
            if problem.branch_e_idx != -1:
                update_eta_sigma(pseudo_plus, pseudo_min, problem.upward_branch, problem.branch_e_idx, problem.parent_lb, problem.parent_branch_e_val, lb)
            # print(f"New suproblem: lb={lb}, ub={ub} (glb={global_lb}, gub={global_ub})")
        except InfeasibleRelaxation:
            # print("Infeasible suproblem...")
            # Pruned by infeasibility.
            continue

        # Update global upper bound.
        if ub < global_ub:
            global_ub = ub
            best_solution = tour
            print(f"Improved upper bound. (glb={global_lb}, gub={global_ub})")

        # by heap property we have that active_nodes[0] has the lowest lower bound
        new_global_lb = min(lb, active_nodes[0].parent_lb) if active_nodes else lb

        # if it's greater than global ub we've already found an optimum and we're just making things worse
        if new_global_lb > global_lb and new_global_lb <= global_ub:
            global_lb = new_global_lb
            print(f"Improved lower bound. (lb={lb} glb={global_lb}, gub={global_ub})")

        # Prune by bound
        if lb > global_ub:
            continue

        # Prune by optimality
        if lb == global_ub:
            continue

        # we use pseudocost branching to find the best variable to branch on
        branch_var_res = find_branch_variable(problem, edges_by_index, pseudo_plus, pseudo_min)
        if branch_var_res is None:
            # we've found an integer solution that the heuristic couldn't find
            print(f"Integer LP solution...")
            if lb < global_ub:
                tour = problem.model.get_tour(edges_by_index)
                global_ub = lb
                best_solution = tour
                print(f"Improved upper bound {global_ub}.")
            continue

        branch_edge, branch_edge_idx, branch_e_value = branch_var_res

        problem_l, problem_r = branch_variable(
            problem, branch_edge, branch_edge_idx, branch_e_value, lb
        )
        heappush(active_nodes, problem_l)
        heappush(active_nodes, problem_r)

    assert global_ub >= global_lb

    opt_time = (perf_counter_ns() - start_opt) / 1e9
    print(f"\t- integer optimal: {length_tour(inst, best_solution)}")
    print(f"\t- optimal tour: {best_solution}")
    timer_times = f"({timer[0]/1e9} s solving LP's. {timer[1]/1e9} s computing heuristics.)"
    print(f"Optimizing using B&B with cutting plane relaxation and Lin-Kernighan took {opt_time} s. {timer_times}")

    # Return optimal solution.
    return global_ub, best_solution


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


def main():
    # inst_path = get_inst_path()
    inst_path = Path("tsp/gr96.dat")
    graph_l = parse_as_adj_matrix(inst_path)
    
    do_branch_and_bound(graph_l)

    gurobi_integer(graph_l)


main()
