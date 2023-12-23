from typing import Literal
from tisp.graph import copy_costs
from tisp.heuristic.lin_kernighan import shuffle_iter, improve_tour
from tisp.types import HeurConstants, Tour, Costs, Edge


def heur_fix(heur_costs: Costs, edge: Edge, fix_val: Literal[0, 1]):
    new_heur_costs = copy_costs(heur_costs)
    if fix_val == 1:
        new_heur_costs[edge[0]][edge[1]] = -float("inf")
        new_heur_costs[edge[1]][edge[0]] = -float("inf")
    elif fix_val == 0:
        new_heur_costs[edge[0]][edge[1]] = float("inf")
        new_heur_costs[edge[1]][edge[0]] = float("inf")
    return new_heur_costs


def lkh(initial_tour: Tour, c: HeurConstants):
    last_tour = initial_tour
    new_tour = initial_tour

    while new_tour is not None:
        untried_tbase = shuffle_iter(new_tour)

        last_tour = new_tour
        new_tour = None

        while len(untried_tbase) > 0:
            tbase = untried_tbase.pop()
            improve_result = improve_tour(
                last_tour, tbase, set(), c
            )

            if improve_result is None:
                continue
            else:
                improved_tour, _ = improve_result

                new_tour = improved_tour
                break

    return last_tour