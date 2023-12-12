import random
from typing import Iterable


Edge = tuple[int, int]
VertTourNghbs = dict[int, tuple[int, int]]
VertNghbs = dict[int, set[int]]

UsedEdges = dict[int, set[Edge]]
SelectedEdges = dict[int, Edge]

def dict_tour(tour: list[int]) -> dict[int, tuple[int, int]]:
    # adjacent edges in the tour
    d: dict[int, tuple[int, int]] = dict()
    final_i = len(tour)-1
    for i in range(len(tour)):
        if i == 0:
            other_i1 = 1
            other_i2 = final_i
        elif i == final_i:
            other_i1 = 0
            other_i2 = final_i-1
        else:
            other_i1 = i+1
            other_i2 = i-1
        
        d[tour[i]] = (tour[other_i1], tour[other_i2])
    
    return d

def random_tour(n: int) -> list[int]:
    return shuffle_iter(range(n))


def shuffle_iter(itrble: Iterable[int]) -> list[int]:
    """Returns a shuffled copy of the iterable as a list."""
    randoms = [(random.random(), i) for i in itrble]
    randoms.sort(key=lambda l: l[0])
    shuffled = list(map(lambda r: r[1], randoms))
    return shuffled

def tour_cost(costs: list[list[float]], tour: list[int]):
    cost = costs[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        cost += costs[tour[i]][tour[i + 1]]

    return cost

def lin_kernighan(costs: list[list[float]]):
    best_cost = float('inf')
    best_tour = None
    runs = 1
    for i in range(runs):
        cost, new_tour = find_tour(costs)
        if cost < best_cost:
            best_cost = cost
            best_tour = new_tour
    
    return best_cost, best_tour

def find_tour(costs: list[list[float]]) -> tuple[float, list[int]]:
    n = len(costs)
    tour =  random_tour(n)

    return improve_tour(tour, costs)


def moves_for_t2(t1: int, t2: int, d_tour: VertTourNghbs, costs: list[list[float]]):
    t_gain = costs[t1][t2]
    while True:
        t2, move_gain, tour_after_move = best_move(t1, t2, t_gain)
        if t2 is None:
            break
        if move_gain > 0:
            return move_gain, tour_after_move
        
    return None


def edge_in_tour(tour: list[int], e: Edge):
    prev = tour[-1]
    for nd in tour:
        if prev == e[0] and nd == e[1]:
            return True
        if prev == e[1] and nd == e[0]:
            return True


def improve_tour(tour: list[int], costs: list[list[float]]) -> tuple[float, list[int]]:
    d_tour = dict_tour(tour)
    cost = tour_cost(costs, tour)
    n = len(costs)
    gain = 0
    current_tour = tour.copy()
    current_d_tour = dict_tour(current_tour)
    

    while True:
        if gain < 0:
            break
        
        for t1 in range(n):
            for x in [0, 1]:
                # pred or succ
                t2 = d_tour[t1][x]

                if edge_in_tour(current_tour, (t1, t2)):
                    continue

                move_res = moves_for_t2(t1, t2, current_d_tour, costs)

                if move_res is None:
                    continue
                else:
                    move_gain, tour_after_move = move_res
                    cost =- move_gain
                    current_tour = tour_after_move
                    current_d_tour = dict_tour(current_tour)
                    break
            
        
        move_gain, tour_after_move = move_23()
        if move_gain > 0:
            cost =- move_gain
            current_tour = tour_after_move
            current_d_tour = dict_tour(current_tour)


    return cost, current_tour


def best_move(t1, t2):
    pass

def move_23() -> :
    pass