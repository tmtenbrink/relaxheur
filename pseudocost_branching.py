import math
import gurobipy as gp


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


PseudoList = tuple[list[float], list[float], list[int]]


def pseudocost(
    x_values: list[float],
    pseudo_plus_list: PseudoList,
    pseudo_min_list: PseudoList,
    mu=1 / 6,
    epsilon=0.0001,
) -> int:
    sigma_j_plus, eta_j_plus, initialized_j_plus = pseudo_plus_list
    sigma_j_min, eta_j_min, initialized_j_min = pseudo_min_list

    best_j = -1
    best_j_score = float("-inf")
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

    return best_j


def update_eta_sigma(
    pseudo_plus_list: PseudoList,
    pseudo_min_list: PseudoList,
    is_upward: bool,
    var_j: int,
    parent_obj: float,
    branch_var: gp.Var,
    obj_val: float,
):
    sigma_j_plus, eta_j_plus, initialized_j_plus = pseudo_plus_list
    sigma_j_min, eta_j_min, initialized_j_min = pseudo_min_list

    if is_upward:
        f_j_plus = math.ceil(branch_var.X) - branch_var.X
        p_i_plus = (obj_val - parent_obj) / f_j_plus

        if eta_j_plus[var_j] == 0:
            initialized_j_plus.append(var_j)
        eta_j_plus[var_j] += 1
        sigma_j_plus[var_j] += p_i_plus
    else:
        f_j_min = branch_var.X - math.floor(branch_var.X)
        p_i_min = (obj_val - parent_obj) / f_j_min

        if eta_j_min[var_j] == 0:
            initialized_j_min.append(var_j)
        eta_j_min[var_j] += 1
        sigma_j_min[var_j] += p_i_min

    return None
