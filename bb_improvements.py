

from branch_bound_all import EdgeValue


def sw_minimum_cut_phase(
    graph: list[list[EdgeValue]], a: int
) -> tuple[int, int, float, list[tuple[int, int]]]:
    graph_n = len(graph)
    A = [a]
    cut_edges: list[tuple[int, int]] = []
    max_cut_weight: float = -1

    while len(A) < graph_n:
        max_cut_weight = -1
        u = 0

        for v in range(graph_n):
            if v not in A:
                cut_weight = sum(graph[v][w][0] for w in A)
                if cut_weight > max_cut_weight:
                    max_cut_weight = cut_weight
                    u = v
                    cut_edges = []
                    for w in A:
                        cut_edges += graph[v][w][1]

        A.append(u)

    s = min(A[-1], A[-2])
    t = max(A[-1], A[-2])

    return s, t, max_cut_weight, cut_edges


def sw_minimum_cut(graph: list[list[EdgeValue]]):
    """Find the minimum cut of a graph using the Stoer-Wagner algorithm."""
    n = len(graph)

    min_cut = float("inf")
    contractions = []
    phase = 0
    best_edge_list = []

    while n > 1:
        a = 0  # Any vertex from V
        s, t, cut_weight, cut_edges = sw_minimum_cut_phase(graph[:n][:n], a)
        if cut_weight < min_cut:
            min_cut = cut_weight
            best_edge_list = cut_edges

        # Merge vertices s and t
        contractions.append((s, t))
        for i in range(n):
            if i != t:
                graph[s][i] = (
                    graph[s][i][0] + graph[t][i][0],
                    graph[s][i][1] + graph[t][i][1],
                )
                graph[i][s] = graph[s][i]

        for i in range(t, n - 1):
            for j in range(n):
                graph[i][j] = graph[i + 1][j]
                graph[j][i] = graph[i][j]

        n -= 1
        phase += 1

    return min_cut, best_edge_list