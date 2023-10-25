import numpy as np
import networkx as nx
import random


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
    stack = [random.randint(0, n)]  # [0]  # Start from vertex 0 -> random

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
    adjacency_matrix = np.array([[0] * n for _ in range(n)])

    for i in range(n):
        for j in range(i + 1, n):
            idx = edge_idx(i, j, n)
            if x_values[idx] == 1:
                adjacency_matrix[i][j] = 1
                adjacency_matrix[j][i] = 1

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


def compute_min_cut(graph: np.ndarray, n: int):
    """Compute min cut using Stoer–Wagner algorithm using Networkx"""
    G = nx.from_numpy_array(graph)

    if not nx.is_connected(G):
        return 0
    # if not is_connected(graph, n):
    #    return 0

    cut_value, partition = nx.stoer_wagner(G)
    return cut_value, partition


n = 6
x_values = np.array([1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1])

# adjacency_matrix = adjacency_matrix_from_vector(x_values, n)
# compute_min_cut(adjacency_matrix, n)

print(check_subtour(x_values, n))
