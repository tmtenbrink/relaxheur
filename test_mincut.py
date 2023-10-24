import numpy as np
import networkx as nx


def edge_idx(lower_i, higher_j, n: int):
    """lower_i and higher_j must be of same dimension or one of them must be a scalar."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)


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
    """Compute min cut using Stoer–Wagner algorithm"""
    if not is_connected(graph, n):
        return 0

    G = nx.from_numpy_array(graph)
    cut_value, partition = nx.stoer_wagner(G)

    return cut_value, partition


def compute_min_cut2(graph: np.ndarray, n: int):
    """Compute min cut using Stoer–Wagner algorithm using Networkx"""
    G = nx.from_numpy_array(graph)
    if not nx.is_connected(G):
        return 0

    cut_value, partition = nx.stoer_wagner(G)
    return cut_value, partition


n = 5
m_edges = (n * (n - 1)) // 2
x_values = np.array([1, 1, 1, 0, 1, 0, 1, 1, 1, 1])
adjacency_matrix = adjacency_matrix_from_vector(x_values, n)

print(is_connected(adjacency_matrix, n))
print(compute_min_cut(adjacency_matrix, n))

print(compute_min_cut2(adjacency_matrix, n))
