from tisp.types import Costs, EdgeValues, EdgesByIndex, Tour


EPSILON = 0.0000001


def normalize_tour(tour: list[int]) -> Tour:
    i_zero = tour.index(0)
    prev_part = tour[:i_zero]
    normalized_tour = tour[i_zero + 1 :] + prev_part
    if normalized_tour[0] < normalized_tour[-1]:
        return [0] + normalized_tour
    else:
        return [0] + list(reversed(normalized_tour))


def tour_from_edge_values(edge_values: EdgeValues, edges_by_index: EdgesByIndex):
    """This fails if the edges are not really a valid tour!"""

    tour = []
    queue = []
    for e_i in range(len(edge_values)):
        x_val = edge_values[e_i]
        if abs(x_val) > EPSILON:
            node_1, node_2 = edges_by_index[e_i]
            queue.append((node_1, node_2))

    while queue:
        # 2 cases
        # - node_1 is already in our tour
        #   in that case, since each edge is only once in our list, either:
        #   + node_1 is at the start
        #   + node_1 is at the end
        #   (it cannot be in the middle because then it already has 2 edges connected and it's an invalid tour)
        # - node_1 is not in the tour
        #   in that case, we do same test for node_2
        #   + node_2 is in the tour
        #     * node_2 is at the start
        #     * node_2 is at the end
        #   + node_2 is not in the tour
        #     we put this edge at the end of the queue and continue until we found one that is there
        node_1, node_2 = queue.pop(0)

        if len(tour) == 0:
            tour += [node_1, node_2]
        elif tour[0] == node_1:
            tour.insert(0, node_2)
        elif tour[-1] == node_1:
            tour.append(node_2)
        elif tour[0] == node_2:
            tour.insert(0, node_1)
        elif tour[-1] == node_2:
            tour.append(node_1)
        else:
            queue.append([node_1, node_2])

    # the start/endpoint is going to be added twice, so we have to remove it
    return normalize_tour(tour[1:])


def tour_cost(graph: Costs, tour: Tour):
    cost = graph[tour[-1]][tour[0]]
    for i in range(len(tour) - 1):
        cost += graph[tour[i]][tour[i + 1]]

    return cost
