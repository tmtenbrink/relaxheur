import unittest
from enum import Enum, unique
import math
from itertools import product
import bisect

from typing import List, Tuple, Optional, Literal

import gurobipy as gp


class UnknownVariableError(Exception):
    pass


# Returns the Models for the left and right branches
def branch_variable(var: gp.Var, m: gp.Model) -> tuple[gp.Model, gp.Model]:
    nm = var.VarName

    left_m = m.copy()
    left_var = left_m.getVarByName(nm)  # Find variable in model
    if left_var is None:
        raise UnknownVariableError(f"Cannot find variable {nm}")
    left_m.addConstr(left_var == 0)  # Add constraint to model

    right_m = m.copy()
    right_var = right_m.getVarByName(nm)  # Find variable in model
    if right_var is None:
        raise UnknownVariableError(f"Cannot find variable {nm}")
    right_m.addConstr(right_var == 1)  # Add constraint to model

    # left is downward, right is upward
    return left_m, right_m


def check_integral(variables: list[gp.Var], epsilon=0.001):
    var_arr = np.array([v.x for v in variables])
    dist = np.abs(np.round(var_arr) - var_arr)
    return np.all(dist < epsilon)


def closest_fractional(variables: list[mip.Var]) -> tuple[mip.Var, int]:
    var_arr = np.array([v.x for v in variables])
    frac, _ = np.modf(var_arr)
    dist_to_05 = np.abs(frac - 0.5)
    closest_to_05_i = np.argmin(dist_to_05, keepdims=True)[0]
    return variables[closest_to_05_i], closest_to_05_i


def first_non_integer(variables: list[gp.Var], epsilon=0.001) -> tuple[gp.Var, int]:
    i = 0
    for v in variables:
        frac = abs(round(v.X) - v.X)
        if frac > epsilon:
            return v, i
        i += 1
    # Return at least some variable
    return variables[0], i


def round_method(x):
    return math.ceil(x)


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


class Solution:
    no_solution: Tuple[List[int], float]
    sigma_j_plus: list[float]
    eta_j_plus: list[float]
    # initialized indexes
    initialized_j_plus: list
    sigma_j_min: list[float]
    eta_j_min: list[float]
    # initialized indexes
    initialized_j_min: list

    def __init__(self):
        self.no_solution = ([], float("inf"))
        self.initialized_j_plus = []
        self.initialized_j_min = []

    def pseudocost(
        self, x_values: list[gp.Var], mu=1 / 6, epsilon=0.0001
    ) -> tuple[gp.Var, int]:
        best_j = -1
        best_j_score = float("-inf")
        for j, v in enumerate(x_values):
            # We definitely do not branch on integer variables
            frac = abs(v.X - round(v.X))
            if frac < epsilon:
                continue

            pseudo_plus = compute_pseudoscore(
                j, self.sigma_j_plus, self.eta_j_plus, self.initialized_j_plus
            )
            pseudo_min = compute_pseudoscore(
                j, self.sigma_j_min, self.eta_j_min, self.initialized_j_min
            )

            f_j_plus = math.ceil(v.X) - v.X
            f_j_min = v.X - math.floor(v.X)

            min_arg = f_j_min * pseudo_min
            plus_arg = f_j_plus * pseudo_plus
            score = (1 - mu) * min(plus_arg, min_arg) + mu * max(plus_arg, min_arg)
            if score > best_j_score:
                best_j = j
                best_j_score = score

        return x_values[best_j], best_j

    def branch_and_bound(self, m: gp.Model) -> Tuple[List[int], float]:
        m.verbose = 0
        m.optimize()

        # Prune by infeasibility
        if m.Status == gp.GRB.INFEASIBLE:
            return self.no_solution

        baseline_value = float(m.ObjVal)

        C = [(baseline_value, m)]

        var_len = len(m.vars)
        self.sigma_j_plus = [0.0] * var_len
        self.eta_j_plus = [0.0] * var_len
        self.sigma_j_min = [0.0] * var_len
        self.eta_j_min = [0.0] * var_len

        upper = float("inf")
        lower = float("-inf")
        best_solution = None
        bound_found = False
        order = -1

        while len(C) > 0:
            c: tuple[float, gp.Model] = C.pop()
            problem = c[1]

            problem_obj = float(problem.objective_value)
            bound_condition_upper = problem_obj < upper
            bound_condition_lower = math.ceil(problem_obj) >= upper

            # Prune by optimality
            if check_integral(problem.vars) and bound_condition_upper:
                if not bound_found:
                    C.sort(key=lambda x: order * x[0])
                    bound_found = True
                upper = problem_obj
                lower = problem_obj

                best_solution = problem
                continue

            # Prune by bound
            elif bound_condition_lower:
                continue

            # Branch the problem
            else:
                # Select which problem/variable to branch on
                branch_var, var_j = self.pseudocost(x_values)

                left_m, right_m = branch_variable(branch_var, problem)

                # Best first search
                left_m.optimize()
                right_m.optimize()

                if not right_m.Status == gp.GRB.INFEASIBLE:
                    right_obj = float(right_m.objective_value)
                    self.update_eta_sigma(
                        True, var_j, problem_obj, branch_var, right_obj
                    )
                    to_append = (right_obj, right_m, True, var_j, problem_obj)
                    if bound_found:
                        bisect.insort(
                            C, (right_obj, right_m), key=lambda x: order * x[0]
                        )
                    else:
                        C.append(to_append)
                if not left_m.status == gp.GRB.INFEASIBLE:
                    left_obj = float(left_m.objective_value)
                    self.update_eta_sigma(
                        False, var_j, problem_obj, branch_var, left_obj
                    )
                    to_append = (left_obj, left_m, False, var_j, problem_obj)
                    if bound_found:
                        bisect.insort(C, (left_obj, left_m), key=lambda x: order * x[0])
                    else:
                        C.append(to_append)

        if best_solution is None:
            return self.no_solution

        return ([x.x for x in best_solution.vars], best_solution.objective_value)

    def update_eta_sigma(
        self,
        is_upward: bool,
        var_j: int,
        parent_obj: float,
        branch_var: gp.Var,
        obj_val: float,
    ):
        if is_upward:
            f_j_plus = math.ceil(branch_var.X) - branch_var.X
            p_i_plus = (obj_val - parent_obj) / f_j_plus

            if self.eta_j_plus[var_j] == 0:
                self.initialized_j_plus.append(var_j)
            self.eta_j_plus[var_j] += 1
            self.sigma_j_plus[var_j] += p_i_plus
        else:
            f_j_min = branch_var.X - math.floor(branch_var.X)
            p_i_min = (obj_val - parent_obj) / f_j_min

            if self.eta_j_min[var_j] == 0:
                self.initialized_j_min.append(var_j)
            self.eta_j_min[var_j] += 1
            self.sigma_j_min[var_j] += p_i_min

        return None
