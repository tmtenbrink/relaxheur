import itertools

import numpy as np
import networkx as nx
import random
from itertools import islice


def edge_idx(lower_i, higher_j, n: int):
    """lower_i and higher_j must be of same dimension or one of them must be a scalar."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)


def get_lt_and_gt_cut(vertex_i: int, n: int):
    # index is based on the lower index vertex
    lt_idx = np.arange(0, vertex_i, dtype=int)
    gt_idx = np.arange(vertex_i + 1, n, dtype=int)

    less_idx = edge_idx(lt_idx, vertex_i, n)
    more_idx = edge_idx(vertex_i, gt_idx, n)

    return less_idx, more_idx


def get_cut_idx(vertex_i: int, n: int):
    less_idx, more_idx = get_lt_and_gt_cut(vertex_i, n)

    cut_idx = np.zeros(n - 1, dtype=int)
    cut_idx[:vertex_i] = less_idx
    cut_idx[vertex_i:] = more_idx

    return cut_idx


def check_subtour(x_values: np.ndarray, n: int):
    """Check if graph, given the char. vector, is connected.
    Start at a random node, then find all nodes that are connected to it.
    Output: #nodes in U and partition U, which are the visited nodes from the initial node.

    Complexity: O(|V|+|E|)"""
    nodes = np.arange(n)
    visited = set()
    stack = [random.randint(0, n - 1)]  # [0]  # Start from vertex 0 -> random

    while stack:
        node = stack.pop()
        visited.add(node)
        cut_idx = get_cut_idx(node, n)
        neighbors = np.delete(nodes, node)
        for i in range(len(cut_idx)):
            if x_values[cut_idx[i]] == 1 and neighbors[i] not in visited:
                stack.append(neighbors[i])

    return len(visited), visited


def adjacency_matrix_from_vector(x_values: np.ndarray, n: int):
    """Create an adjacency matrix from x_values"""
    if len(x_values) != (n * (n - 1)) // 2:
        raise ValueError("x_values does not have the correct length")

    adjacency_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            idx = edge_idx(i, j, n)
            adjacency_matrix[i][j] = x_values[idx]
            adjacency_matrix[j][i] = adjacency_matrix[i][j]

    return adjacency_matrix


def is_connected(graph: np.ndarray, n: int):
    """Manually check if graph, given its adjaceny matrix, is connected.
    Complexity: O(|V|+|E|)"""
    visited = set()
    stack = [0]  # Start from vertex 0

    while stack:
        node = stack.pop()
        visited.add(node)
        for neighbor in range(n):
            if graph[node][neighbor] == 1 and neighbor not in visited:
                stack.append(neighbor)

    return len(visited) == n


def get_edge_tpls_arr(n: int):
    tpls = []
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            tpls.append((i, j, 0))

            index += 1

    return np.array(tpls)


def arr_to_edge_tpls(arr: np.ndarray):
    return list(map(tuple, arr))


def partition_to_edge_idx(n: int, partition: tuple[list, list]):
    # unfortunately networkx only returns a partition, not the edges in the cut
    vertices_one = np.array(partition[0])
    vertices_two = np.array(partition[1])
    # we want all possible combinations, as that is the edges of our cut as it is fully connected
    # meshgrid will repeat first array in along columns, second along rows, returning two arrays
    # so if you then take any element and pair it with the same location in the second array you get a combination
    grid = np.meshgrid(vertices_one, vertices_two)
    # we then stack them along the third dimension as described, so we the values at each point are a combiantion
    # of the two arrays
    stacked = np.dstack(grid)
    # finally we reshape the array into a simple 2D array, where each row is an edge
    # the -1 indicates that the shape of the array in the first dimension is inferred (the number of edges)
    edges = stacked.reshape(-1, 2)
    # we sort the edges so that the lowest index is first
    edges = np.sort(edges, axis=1)
    # finally we compute the edge index of each edge
    return edge_idx(edges[:, 0], edges[:, 1], n)


def compute_min_cut(
    x_values: np.ndarray, n: int, base_graph: nx.Graph, edge_tuples_arr: np.ndarray
):
    """Compute min cut using Stoer–Wagner algorithm using Networkx. x_values should be 1D array with our edge
    index convention. edge_tuples_arr should be of the same length, but here each element is a 3-element array
    (so an m*3 array) that is (i, j, w) with i the lower vertex index, j the higher and w to be used for weight.
    """
    # We make a copy so we do not modify the base graph
    weighted_graph = base_graph.copy()
    # We copy the array with weights zero
    edge_tuples_arr = edge_tuples_arr.copy()
    edge_tuples_arr[:, 2] = x_values
    weighted_graph.add_weighted_edges_from(edge_tuples_arr)

    cut_value, partition = nx.stoer_wagner(weighted_graph)
    cut_edges = partition_to_edge_idx(n, partition)

    return cut_value, cut_edges


def get_cut_idx_arr(vertices: np.ndarray, n: int):
    """Vertices should be 1D array."""
    vertices_len = vertices.size
    # each column is the vertex repeated
    vertices_repeated = np.tile(vertices, (n - 1, 1))
    # each column is the index from 0 to n-1
    indexes_repeated = np.repeat(np.arange(n - 1), vertices_len).reshape(
        (n - 1, vertices_len)
    )
    # boolean mask of elements where index is less than vertex
    lt_mask = indexes_repeated < vertices_repeated
    # boolean mask of where index is greater/equal than vertex
    gt_mask = indexes_repeated >= vertices_repeated
    # we add one to the ones that are greater or equal, so that vertex itself doesn't show up
    # now each column is all the vertices not equal to the vertex of that column
    indexes_repeated[gt_mask] += 1

    # we can now use edge_idx as before
    vertices_repeated[gt_mask] = edge_idx(
        vertices_repeated[gt_mask], indexes_repeated[gt_mask], n
    )
    vertices_repeated[lt_mask] = edge_idx(
        indexes_repeated[lt_mask], vertices_repeated[lt_mask], n
    )

    return vertices_repeated


def cut_out(vertices: np.ndarray, n: int):
    all_edges = get_cut_idx_arr(vertices, n)

    return all_edges * 2


def cut_in(vertices: np.ndarray, n: int):
    all_edges = get_cut_idx_arr(vertices, n)

    return all_edges * 2 + 1


def sw_minimum_cut_phase(graph, a):
    n = len(graph)
    A = [a]
    w_list = []

    while len(A) < n:
        max_cut_weight = -1

        for v in range(n):
            if v not in A:
                cut_weight = sum(graph[v][w][0] for w in A)
                if cut_weight > max_cut_weight:
                    max_cut_weight = cut_weight
                    u = v
                    w_list = []
                    for w in A:
                        w_list += graph[v][w][1]

        A.append(u)

    s = min(A[-1], A[-2])
    t = max(A[-1], A[-2])

    return s, t, cut_weight, w_list


def sw_minimum_cut(graph):
    n = len(graph)

    min_cut = float("inf")
    contractions = []
    phase = 0
    best_w_list = []

    while n > 1:
        a = 0  # Any vertex from V
        s, t, cut_weight, w_list = sw_minimum_cut_phase(graph[:n][:n], a)
        if cut_weight < min_cut:
            min_cut = cut_weight
            best_phase = phase
            best_w_list = w_list

        # Merge vertices s and t
        contractions.append((s, t))
        print(s, t, cut_weight)
        for i in range(n):
            if i != t:
                graph[s][i] = [
                    graph[s][i][0] + graph[t][i][0],
                    graph[s][i][1] + graph[t][i][1],
                ]
                graph[i][s] = graph[s][i]

        for i in range(t, n - 1):
            for j in range(n):
                graph[i][j] = graph[i + 1][j]
                graph[j][i] = graph[i][j]

        n -= 1
        phase += 1
    # Recover partition
    print(contractions)
    for i in islice(contractions, best_phase):
        print(i)

    v = contractions[best_phase][1]

    return min_cut, phase, best_w_list


graph = [
    [
        [0, [(0, 0)]],
        [5, [(0, 1)]],
        [0, [(0, 2)]],
        [0, [(0, 3)]],
        [1, [(0, 4)]],
        [4, [(0, 5)]],
    ],
    [
        [5, [(1, 0)]],
        [0, [(1, 1)]],
        [2, [(1, 2)]],
        [0, [(1, 3)]],
        [0, [(1, 4)]],
        [0, [(1, 5)]],
    ],
    [
        [0, [(2, 0)]],
        [2, [(2, 1)]],
        [0, [(2, 2)]],
        [6, [(2, 3)]],
        [1, [(2, 4)]],
        [1, [(2, 5)]],
    ],
    [
        [0, [(3, 0)]],
        [0, [(3, 1)]],
        [6, [(3, 2)]],
        [0, [(3, 3)]],
        [3, [(3, 4)]],
        [0, [(3, 5)]],
    ],
    [
        [1, [(4, 0)]],
        [0, [(4, 1)]],
        [1, [(4, 2)]],
        [3, [(4, 3)]],
        [0, [(4, 4)]],
        [0, [(4, 5)]],
    ],
    [
        [4, [(5, 0)]],
        [0, [(5, 1)]],
        [1, [(5, 2)]],
        [0, [(5, 3)]],
        [0, [(5, 4)]],
        [0, [(5, 5)]],
    ],
]
# G = nx.from_numpy_array(graph)
# print(nx.is_connected(G))
# cut_value, partition = nx.stoer_wagner(G)
# print(cut_value, partition)

print(sw_minimum_cut(graph))
