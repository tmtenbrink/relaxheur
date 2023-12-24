from heapq import heapify, heappop, heappush
from time import perf_counter_ns
from tisp.branch_bound.branching import (
    branch_variable,
    find_branch_variable,
    update_eta_sigma,
)
from tisp.error import InfeasibleRelaxation
from tisp.graph import best_neighbors, copy_costs
from tisp.heuristic.compute import lkh
from tisp.lp.model import lp_edges, lp_initialize
from tisp.lp.optimize import lp_optimize_optimal
from tisp.tour import tour_cost, tour_from_edge_values

from tisp.types import (
    Costs,
    Edge,
    Formulation,
    HeurCosts,
    LPConstants,
    LPModel,
    Subproblem,
    SubproblemState,
    Timer,
)


def compute_bounds(
    costs: list[list[float]],
    problem: SubproblemState,
    best_solution: list[int],
    best_solution_cost: float,
    timer: Timer,
):
    lp, _, _, _, _, _, heur_cost = problem

    # Solve relaxation to find LB
    # RAISES InfeasibleRelaxation
    before = perf_counter_ns()
    lb = lp_optimize_optimal(lp)
    after_lp = perf_counter_ns()

    heur_tour, heur_ub = heuristic(costs, heur_cost, best_solution, best_solution_cost)

    after_heur = perf_counter_ns()

    timer = (timer[0] + (after_lp - before), timer[1] + (after_heur - after_lp))

    return lb, heur_ub, heur_tour, timer


def bb_subproblem_pack(
    lp: LPModel,
    fixed_one: list[Edge],
    fixed_zero: list[Edge],
    branch_edge: int,
    is_upward: bool,
    branch_edge_val: float,
    heur_costs: HeurCosts,
) -> SubproblemState:
    return (
        lp,
        fixed_one,
        fixed_zero,
        branch_edge,
        is_upward,
        branch_edge_val,
        heur_costs,
    )


def bb_intialize(
    costs: Costs, formulation=Formulation.CUTTING_PLANE
) -> tuple[Subproblem, LPConstants]:
    lp = lp_initialize(formulation, costs)
    _, _, c = lp

    return Subproblem(
        -1, bb_subproblem_pack(lp, [], [], -1, False, -1, copy_costs(costs))
    ), c


def heuristic(
    costs: Costs,
    heur_cost: HeurCosts,
    best_solution: list[int],
    best_sol_cost: float,
) -> tuple[list[int], float]:
    n = len(costs)
    best_nbs = best_neighbors(heur_cost)

    lkh_tour = lkh(best_solution, (n, heur_cost, best_nbs))
    lkh_cost = tour_cost(costs, lkh_tour)

    return lkh_tour, lkh_cost


def initial_ub(inst: Costs) -> float:
    ub = 0.0
    # cost is definitely lower than the sum of all the costs
    for lst in inst:
        ub += sum(lst)
    return ub


def do_branch_and_bound(inst: Costs):
    global_problem, c = bb_intialize(inst)
    n, m_edges, _, _, edges_by_index = c

    global_lb = global_problem.parent_lb

    base_tour = list(range(n))
    base_cost = tour_cost(inst, base_tour)

    best_solution, global_ub = heuristic(
        inst, global_problem.state[6], base_tour, base_cost
    )

    pseudo_plus = ([0.0] * m_edges, [0.0] * m_edges, list[int]())
    pseudo_min = ([0.0] * m_edges, [0.0] * m_edges, list[int]())

    # Initialization.
    active_nodes: list[Subproblem] = [global_problem]
    # we use a priority queue since we use the heuristic to quickly get good upper bounds
    # this is best-first search
    heapify(active_nodes)

    start_opt = perf_counter_ns()
    timer: Timer = (0, 0)

    while active_nodes:
        # Select an active node to process. We choose the one with the lowest lower bound (hopefully this will also
        # have a good upper bound, alllowing us to prune faster)
        problem = heappop(active_nodes)

        lp, _, _, branch_edge, is_upward, branch_edge_value, _ = problem.state

        try:
            lb, ub, heur_tour, timer = compute_bounds(
                inst, problem.state, best_solution, global_ub, timer
            )

            if branch_edge != -1:
                update_eta_sigma(
                    pseudo_plus,
                    pseudo_min,
                    is_upward,
                    branch_edge,
                    problem.parent_lb,
                    branch_edge_value,
                    lb,
                )
        except InfeasibleRelaxation:
            # Pruned by infeasibility.
            continue

        # Update global upper bound.
        if ub < global_ub:
            global_ub = ub
            best_solution = heur_tour
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

        edge_values, _ = lp_edges(lp)
        # we use pseudocost branching to find the best variable to branch on
        branch_var_res = find_branch_variable(
            edge_values, edges_by_index, pseudo_plus, pseudo_min
        )
        if branch_var_res is None:
            # we've found an integer solution that the heuristic couldn't find
            print(f"Integer LP solution...")
            if lb < global_ub:
                tour = tour_from_edge_values(edge_values, edges_by_index)
                global_ub = lb
                best_solution = tour
                print(f"Improved upper bound {global_ub}.")
            continue

        branch_edge, branch_edge_idx, branch_e_value = branch_var_res

        problem_l, problem_r = branch_variable(
            problem.state, branch_edge, branch_edge_idx, branch_e_value, lb
        )
        heappush(active_nodes, problem_l)
        heappush(active_nodes, problem_r)

    assert global_ub >= global_lb

    opt_time = (perf_counter_ns() - start_opt) / 1e9
    print(f"\t- integer optimal: {tour_cost(inst, best_solution)}")
    print(f"\t- optimal tour: {best_solution}")
    timer_times = (
        f"({timer[0]/1e9} s solving LP's. {timer[1]/1e9} s computing heuristics.)"
    )
    print(
        f"Optimizing using B&B with cutting plane relaxation and Lin-Kernighan took {opt_time} s. {timer_times}"
    )

    # Return optimal solution.
    return global_ub, best_solution
