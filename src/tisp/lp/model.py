from numbers import Real
from time import perf_counter_ns
from typing import Literal, cast
import mip

from tisp.error import UndefinedVariableError
from tisp.graph import edge_idxs_for_all_v, get_edges_by_index
from tisp.types import (
    Costs,
    EdgeCosts,
    EdgeValues,
    EdgesByIndex,
    Formulation,
    LPConstants,
    LPModel,
    VertexEdges,
)


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


def lp_pack_const(
    n: int,
    m_edges: int,
    vertex_edges: VertexEdges,
    edge_costs: EdgeCosts,
    edges_by_index: EdgesByIndex,
) -> LPConstants:
    return n, m_edges, vertex_edges, edge_costs, edges_by_index


def cutting_plane_model(m_edges: int, edge_costs: EdgeCosts):
    model = mip.Model("cutting_plane", solver_name=mip.CBC)
    model.verbose = 0
    model.lp_method = mip.LP_Method.DUAL
    model.emphasis = mip.SearchEmphasis.OPTIMALITY

    print("Setting up cutting plane model...")
    start_build = perf_counter_ns()

    z_expr = dict(
        [
            (
                model.add_var(lb=cast(Real, 0.0), ub=cast(Real, 1.0), name=f"edge_{e}"),
                cast(Real, edge_costs[e]),
            )
            for e in range(m_edges)
        ]
    )
    z = mip.LinExpr(expr=z_expr)
    model.objective = mip.minimize(z)

    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return model


def lp_initialize(formulation: Formulation, costs: Costs) -> LPModel:
    edge_costs = compute_edge_costs(costs)
    n = len(costs)
    m_edges = (n * (n - 1)) // 2

    if formulation == Formulation.CUTTING_PLANE:
        model = cutting_plane_model(m_edges, edge_costs)
    else:
        raise NotImplementedError()

    vertex_edges = edge_idxs_for_all_v(n)
    edges_by_index = get_edges_by_index(n)
    c = lp_pack_const(n, m_edges, vertex_edges, edge_costs, edges_by_index)
    return model, formulation, c


def lp_copy(lp: LPModel) -> LPModel:
    model, formulation, c = lp
    copied_model = model.copy()
    copied_model.verbose = model.verbose

    return copied_model, formulation, c


def lp_copy_fix(lp: LPModel, edge_idx: int, fix_val: Literal[0, 1]) -> LPModel:
    copy_model, formulation, c = lp_copy(lp)
    edge_var = copy_model.var_by_name("edge_" + str(edge_idx))
    copy_model += edge_var == fix_val, f"x_{edge_idx}=={fix_val}"

    return copy_model, formulation, c


def lp_edges(lp: LPModel) -> tuple[EdgeValues, list[mip.Var]]:
    model, _, c = lp
    _, m, _, _, _ = c
    edge_values: list[float] = []
    edge_vars: list[mip.Var] = []
    for e in range(m):
        v = model.var_by_name(f"edge_{e}")
        if v is None or (v_x := v.x) is None:
            raise UndefinedVariableError()

        edge_values.append(float(v_x))
        edge_vars.append(v)

    return edge_values, edge_vars
