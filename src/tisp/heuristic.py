from typing import Literal
from tisp.graph import copy_costs

from tisp.types import Costs, Edge


def heur_fix(heur_costs: Costs, edge: Edge, fix_val: Literal[0, 1]):
    new_heur_costs = copy_costs(heur_costs)
    if fix_val == 1:
        new_heur_costs[edge[0]][edge[1]] = -float("inf")
        new_heur_costs[edge[1]][edge[0]] = -float("inf")
    elif fix_val == 0:
        new_heur_costs[edge[0]][edge[1]] = float("inf")
        new_heur_costs[edge[1]][edge[0]] = float("inf")
    return new_heur_costs
