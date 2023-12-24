import math
from typing import Literal, Optional
from tisp.heuristic.compute import heur_fix
from tisp.lp.model import lp_copy_fix

from tisp.types import Edge, EdgeValues, PseudoList, Subproblem, SubproblemState


def compute_pseudoscore(
    index_j: int, sigma_j: list[float], eta_j: list[float], initialized_j: list[int]
):
    uninitialized = eta_j[index_j] == 0
    if uninitialized:
        if not initialized_j:
            return 1

        sigmas_d_etas = [sigma_j[j] / eta_j[j] for j in initialized_j]

        return sum(sigmas_d_etas) / len(initialized_j)

    return sigma_j[index_j] / eta_j[index_j]


def pseudocost(
    x_values: list[float],
    pseudo_plus_list: PseudoList,
    pseudo_min_list: PseudoList,
    mu=1 / 6,
    epsilon=0.0001,
) -> tuple[int, float]:
    sigma_j_plus, eta_j_plus, initialized_j_plus = pseudo_plus_list
    sigma_j_min, eta_j_min, initialized_j_min = pseudo_min_list

    best_j = -1
    best_j_score = float("-inf")
    best_j_value = -1
    for j, value in enumerate(x_values):
        # We definitely do not branch on integer variables
        frac = abs(value - round(value))
        if frac < epsilon:
            continue

        pseudo_plus = compute_pseudoscore(
            j, sigma_j_plus, eta_j_plus, initialized_j_plus
        )
        pseudo_min = compute_pseudoscore(j, sigma_j_min, eta_j_min, initialized_j_min)

        f_j_plus = math.ceil(value) - value
        f_j_min = value - math.floor(value)

        min_arg = f_j_min * pseudo_min
        plus_arg = f_j_plus * pseudo_plus
        score = (1 - mu) * min(plus_arg, min_arg) + mu * max(plus_arg, min_arg)
        if score > best_j_score:
            best_j = j
            best_j_score = score
            best_j_value = value

    return best_j, best_j_value


def update_eta_sigma(
    pseudo_plus_list: PseudoList,
    pseudo_min_list: PseudoList,
    is_upward: bool,
    var_j: int,
    parent_obj: float,
    parent_branch_var_val: float,
    obj_val: float,
):
    sigma_j_plus, eta_j_plus, initialized_j_plus = pseudo_plus_list
    sigma_j_min, eta_j_min, initialized_j_min = pseudo_min_list

    if is_upward:
        f_j_plus = math.ceil(parent_branch_var_val) - parent_branch_var_val
        p_i_plus = (obj_val - parent_obj) / f_j_plus

        if eta_j_plus[var_j] == 0:
            initialized_j_plus.append(var_j)
        eta_j_plus[var_j] += 1
        sigma_j_plus[var_j] += p_i_plus
    else:
        f_j_min = parent_branch_var_val - math.floor(parent_branch_var_val)
        p_i_min = (obj_val - parent_obj) / f_j_min

        if eta_j_min[var_j] == 0:
            initialized_j_min.append(var_j)
        eta_j_min[var_j] += 1
        sigma_j_min[var_j] += p_i_min

    return None


def branch_variable(
    problem: SubproblemState,
    edge: Edge,
    e_idx: int,
    branch_e_val: float,
    parent_lb: float,
):
    lp, fixed_one, fixed_zero, _, _, _, heur_costs = problem

    new_model_l = lp_copy_fix(lp, e_idx, 1)
    new_model_r = lp_copy_fix(lp, e_idx, 0)

    fixed_one_l = fixed_one + [edge]
    fixed_zero_r = fixed_zero + [edge]

    heur_costs_l = heur_fix(heur_costs, edge, 1)
    heur_costs_r = heur_fix(heur_costs, edge, 0)

    problem_l = Subproblem(
        parent_lb,
        (
            new_model_l,
            fixed_one_l,
            fixed_zero,
            e_idx,
            True,
            branch_e_val,
            heur_costs_l,
        ),
    )
    problem_r = Subproblem(
        parent_lb,
        (
            new_model_r,
            fixed_one,
            fixed_zero_r,
            e_idx,
            False,
            branch_e_val,
            heur_costs_r,
        ),
    )

    return problem_l, problem_r


def find_branch_variable(
    edge_values: EdgeValues,
    edges_by_index: dict[int, tuple[int, int]],
    pseudo_plus: PseudoList,
    pseudo_min: PseudoList,
    method: Literal["pseudocost", "first-non-integer"] = "pseudocost",
) -> Optional[tuple[Edge, int, float]]:
    if method == "first-non-integer":
        for e_idx, e_val in enumerate(edge_values):
            epsilon = 0.0000001
            if abs(round(e_val) - e_val) > epsilon:
                # print(f"edge {e_idx} with value {e_val} is not integer. Splitting...")
                return edges_by_index[e_idx], e_idx, e_val
    elif method == "pseudocost":
        e_idx, e_val = pseudocost(edge_values, pseudo_plus, pseudo_min)
        if e_idx == -1:
            return None
        return edges_by_index[e_idx], e_idx, e_val

    return None
