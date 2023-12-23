from time import perf_counter_ns
import mip

from tisp.types import Formulation
from tisp.error import InfeasibleRelaxation, NotOptimalError
from tisp.lp.cut_relax import separation
from tisp.lp.model import LPModel, lp_edges


def cut_lp_optimize(m: mip.Model):
    m.optimize()


def cut_relax_optimize(lp: LPModel) -> mip.Model:
    """RAISES InfeasibleRelaxation when unable to optimize."""
    m, _, c = lp
    
    invalid = True

    while invalid:
        cut_lp_optimize(m)
        if m.status == mip.OptimizationStatus.INFEASIBLE:
            raise InfeasibleRelaxation()
        elif m.status != mip.OptimizationStatus.OPTIMAL:
            raise NotOptimalError("Unexpected non-optimal status!")

        edge_values, edge_vars = lp_edges(lp)
        # this modifies the underlying model and adds constraints
        in_subtour = separation(m, edge_values, edge_vars, c)

        invalid = not in_subtour

    return m


def lp_optimize_optimal(lp: LPModel, show_time=False) -> float:
    start_opt = perf_counter_ns()

    model, formulation, _ = lp
    if formulation == Formulation.CUTTING_PLANE:
        model = cut_relax_optimize(lp)
    else:
        model.optimize()

    
    model_obj = model.objective_value
    if model.status != mip.OptimizationStatus.OPTIMAL or model_obj is None:
        raise NotOptimalError()
    
    opt_time = (perf_counter_ns() - start_opt) / 1e9
    if show_time:
        print(f"Took {opt_time} s.")

    
    return float(model_obj)