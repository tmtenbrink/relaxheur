import random
from typing import Iterable, Optional
from tisp.tour import normalize_tour

from tisp.types import HeurConstants, Edge, Tour


def random_tour(n: int) -> list[int]:
    return shuffle_normalized(list(range(n)))


def shuffle_iter(itrble: Iterable[int]) -> list[int]:
    """Returns a shuffled copy of the iterable as a list."""
    randoms = [(random.random(), i) for i in itrble]
    randoms.sort(key=lambda l_el: l_el[0])
    shuffled = list(map(lambda r: r[1], randoms))
    return shuffled


def shuffle_normalized(tour: Tour) -> list[int]:
    """Returns a shuffled copy of the iterable as a list."""
    shuffled_others = shuffle_iter(tour[1:])
    if shuffled_others[0] < shuffled_others[-1]:
        return [0] + shuffled_others
    else:
        return [0] + list(reversed(shuffled_others))


# def tour_cost(costs: Costs, tour: Tour) -> float:
#     cost = costs[tour[-1]][tour[0]]
#     for i in range(len(tour) - 1):
#         cost += costs[tour[i]][tour[i + 1]]

#     return cost


# def tour_unique(tour: list[int]):
#     max_size = len(tour)

#     bytes_per_n = (max_size.bit_length() + 7) // 8

#     return reduce(lambda b, t: b + t.to_bytes(bytes_per_n, byteorder="big"), tour, b"")


def tour_succ_i(tour: Tour, node_i: int):
    if node_i == len(tour) - 1:
        return tour[0], 0

    return tour[node_i + 1], node_i + 1


def tour_succ(tour: Tour, node: int):
    node_i = tour.index(node)
    return tour_succ_i(tour, node_i)


def tour_pred_i(tour: Tour, node_i: int):
    if node_i == 0:
        return tour[-1], len(tour) - 1

    return tour[node_i - 1], node_i - 1


def tour_pred(tour: Tour, node: int):
    node_i = tour.index(node)
    return tour_pred_i(tour, node_i)


def edge_in_edge_list(e: Edge, edge_list: Iterable[Edge]) -> bool:
    return e in edge_list or (e[1], e[0]) in edge_list


def edge_in_tour(tour: list[int], e: Edge):
    prev = tour[-1]
    for nd in tour:
        if prev == e[0] and nd == e[1]:
            return True
        if prev == e[1] and nd == e[0]:
            return True


def improve_tour(tour: Tour, t0: int, seen_tours: set[bytes], c: HeurConstants):
    _, costs, _ = c
    current_tour = tour.copy()
    t0_i = tour.index(t0)
    for t1, t1_i in (tour_succ_i(tour, t0_i), tour_pred_i(tour, t0_i)):
        broken_cost = costs[t0][t1]

        move_res = best_move(t0_i, t1_i, broken_cost, current_tour, seen_tours, c)

        if move_res:
            return move_res


def swap_any_not_equal(t0_i: int, t1_i: int, t2_i: int, t3_i: int):
    return (
        t0_i != t1_i
        and t0_i != t2_i
        and t0_i != t3_i
        and t1_i != t2_i
        and t1_i != t3_i
        and t2_i != t3_i
    )


def is_swap_feasible(tour: Tour, t0_i: int, t1_i: int, t2_i: int, t3_i: int):
    if not swap_any_not_equal(t0_i, t1_i, t2_i, t3_i):
        return False

    if tour_succ_i(tour, t0_i)[1] == t1_i:
        if t3_i != tour_pred_i(tour, t2_i)[1]:
            return False
    elif tour_pred_i(tour, t0_i)[1] == t1_i:
        if t3_i != tour_succ_i(tour, t2_i)[1]:
            return False
    else:
        raise ValueError("Edge does not exist!")

    return True


def is_swap_feasible_by_v(tour: Tour, t0: int, t1: int, t2: int, t3: int):
    return is_swap_feasible(
        tour, tour.index(t0), tour.index(t1), tour.index(t2), tour.index(t3)
    )


MoveState = tuple[float, int, int, int, int, Tour, set[Edge], set[Edge]]
XResult = Optional[tuple[bool, float, Tour, set[Edge], set[Edge]]]


def select_x(
    gain: float,
    t0_i: int,
    t1_i: int,
    t2_i: int,
    t3_i: int,
    current_tour: Tour,
    broken_edges: set[Edge],
    joined_edges: set[Edge],
    c: HeurConstants,
) -> XResult:
    _, costs, _ = c
    t2 = current_tour[t2_i]
    t3 = current_tour[t3_i]
    broken_edge = (t2, t3)
    broken_cost = costs[t2][t3]

    if t0_i == t3_i:
        return None

    if edge_in_edge_list(broken_edge, broken_edges) or edge_in_edge_list(
        broken_edge, joined_edges
    ):
        return None

    if not is_swap_feasible(current_tour, t0_i, t1_i, t2_i, t3_i):
        return None

    broken_edges.add(broken_edge)

    t0 = current_tour[t0_i]
    joined_edge = (t0, t3)
    joined_cost = costs[t0][t3]
    current_gain = gain + (broken_cost - joined_cost)

    swapped_tour = apply_swap(current_tour, t0_i, t1_i, t2_i, t3_i)

    if current_gain > 0.00001:
        joined_edges.add(joined_edge)
        return True, current_gain, swapped_tour, broken_edges, joined_edges
    else:
        # todo add broken edges and stuff
        return False, current_gain, swapped_tour, broken_edges, joined_edges


YResult = Optional[tuple[float, int, int, int, int, set[Edge], set[Edge]]]


def select_y(
    gain: float,
    t0_i: int,
    t3_i: int,
    current_tour: Tour,
    broken_edges: set[Edge],
    joined_edges: set[Edge],
    c: HeurConstants,
) -> YResult:
    _, costs, best_nbs = c
    t0 = current_tour[t0_i]
    t3 = current_tour[t3_i]

    broken_cost = costs[t0][t3]

    for nd in best_nbs[t3]:
        nd_i = current_tour.index(nd)
        for nd_nb_i, _ in (
            tour_succ_i(current_tour, nd_i),
            tour_pred_i(current_tour, nd_i),
        ):
            if not is_swap_feasible(current_tour, t0_i, t3_i, nd_i, nd_nb_i):
                continue

            joined_edge = (t3, nd)
            joined_cost = costs[t3][nd]

            curr_gain = gain + (broken_cost - joined_cost)

            if curr_gain <= 0.00001:
                continue

            if edge_in_edge_list(joined_edge, joined_edges) or edge_in_edge_list(
                joined_edge, broken_edges
            ):
                continue

            return curr_gain, t0_i, t3_i, nd_i, nd_nb_i, broken_edges, joined_edges

    return None


def apply_swap(tour: Tour, t0_i: int, t1_i: int, t2_i: int, t3_i: int):
    # previously, t0 and t1 are connected; t2 and t3 are connected
    # after the swap, t2 connects to t1 and t3 connects to t0

    sorted_by_i = sorted([t0_i, t1_i, t2_i, t3_i])

    min_i = sorted_by_i[0]
    min_i_1 = sorted_by_i[1]
    min_i_2 = sorted_by_i[2]
    min_i_3 = sorted_by_i[3]

    split_first = (
        len(tour) - 1 == min_i_3
        and min_i == 0
        and (min_i_1 != 1 or (min_i_1 == 1 and min_i_2 == 2))
    )

    if split_first:
        new_tour = (
            [0] + tour[min_i_2 : min_i_3 + 1] + list(reversed(tour[1 : min_i_1 + 1]))
        )
    else:
        new_tour = (
            tour[: min_i + 1]
            + list(reversed(tour[min_i_1 : min_i_2 + 1]))
            + tour[min_i_3:]
        )

    return normalize_tour(new_tour)


class SeenTour(Exception):
    pass


def move_sequence(
    base_start: MoveState, seen_tours: set[bytes], c: HeurConstants
) -> Optional[tuple[Tour, float]]:
    queue: list[MoveState] = [base_start]

    while queue:
        (
            gain,
            t0_i,
            t1_i,
            t2_i,
            t3_i,
            current_tour,
            broken_edges,
            joined_edges,
        ) = queue.pop()

        t0 = current_tour[t0_i]
        t3 = current_tour[t3_i]

        # sw_ind = (t0_i, t1_i, t2_i, t3_i)
        # swap = (t0, current_tour[t1_i], current_tour[t2_i], t3)

        x_res = select_x(
            gain, t0_i, t1_i, t2_i, t3_i, current_tour, broken_edges, joined_edges, c
        )

        if x_res is None:
            continue

        x_improved, gain, new_tour, broken_edges, joined_edges = x_res

        if x_improved:
            return new_tour, gain

        # if tour_unique(new_tour) in seen_tours:
        #     raise SeenTour()

        mod_t0_i = new_tour.index(t0)
        mod_t3_i = new_tour.index(t3)

        y_res = select_y(
            gain, mod_t0_i, mod_t3_i, new_tour, broken_edges, joined_edges, c
        )

        if y_res is None:
            continue

        gain, new_t0_i, new_t1_i, new_t2_i, new_t3_i, broken_edges, joined_edges = y_res

        # new_sw_ind = (t0_i, t1_i, t2_i, t3_i)
        # new_swap = (new_tour[new_t0_i], new_tour[t1_i], new_tour[t2_i], new_tour[t3_i])

        queue.append(
            (
                gain,
                new_t0_i,
                new_t1_i,
                new_t2_i,
                new_t3_i,
                new_tour,
                broken_edges,
                joined_edges,
            )
        )

    return None


def best_move(
    t0_i: int,
    t1_i: int,
    broken_cost: float,
    current_tour: Tour,
    seen_tours: set[bytes],
    constants: HeurConstants,
) -> Optional[tuple[Tour, float]]:
    t1 = current_tour[t1_i]
    _, costs, best_nbs = constants
    for t2 in best_nbs[t1]:
        t2_i = current_tour.index(t2)
        for _, t3_i in (
            tour_succ_i(current_tour, t2_i),
            tour_pred_i(current_tour, t2_i),
        ):
            if not is_swap_feasible(current_tour, t0_i, t1_i, t2_i, t3_i):
                continue

            joined_edge = (t1, t2)
            joined_cost = costs[t1][t2]

            gain = broken_cost - joined_cost

            if gain <= 0.00001:
                continue
            if edge_in_tour(current_tour, joined_edge):
                continue

            t0 = current_tour[t0_i]
            broken_edges = {(t0, t1)}
            joined_edges = {joined_edge}

            maybe_tour_res = move_sequence(
                (
                    gain,
                    t0_i,
                    t1_i,
                    t2_i,
                    t3_i,
                    current_tour,
                    broken_edges,
                    joined_edges,
                ),
                seen_tours,
                constants,
            )

            if maybe_tour_res:
                maybe_tour, gain = maybe_tour_res
                return maybe_tour, gain

    return None
