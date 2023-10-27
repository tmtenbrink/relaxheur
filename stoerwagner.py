import numpy as np
import random
from heapq import heappush, heappop

def edge_idx(lower_i, higher_j, n: int):
    """lower_i and higher_j must be of same dimension or one of them must be a scalar."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)


def search_list(h: list, index: int):
    for i in h:
        if i[0] == index:
            return i[1]
    return 0

def stoer_wagner(x_values: np.ndarray, edge_tuples_arr: np.ndarray, n: int, cut_arr: np.ndarray):
    cut_value = float("inf")
    nodes = set(range(n))
    contractions = []  # contracted node pairs

    # Repeatedly pick a pair of nodes to contract until only one node is left.
    for i in range(n - 1):
        # Pick an arbitrary node u and create a set A = {u}.
        u = random.randint(0, n-1)
        A = {u}
        # Repeatedly pick the node "most tightly connected" to A and add it to
        # A. The tightness of connectivity of a node not in A is defined by the
        # of edges connecting it to nodes in A.
        h = []
        cut = cut_arr[u]
        for e in cut:
            edge = edge_tuples_arr[e]
            other_v = edge[0] if edge[0] != u else edge[1]
            heappush(other_v, (other_v, x_values[e]))
        # for v, e in G[u].items():
        #     h.insert(v, -e["weight"])
        # Repeat until all but one node has been added to A.
        for j in range(n - i - 2):
            u = heappop(h)[0]
            A.add(u)
            cut = cut_arr[u]
            for e in cut:
                edge = edge_tuples_arr[e]
                other_v = edge[0] if edge[0] != u else edge[1]
                if other_v not in A:
                    w = search_list(h, other_v)
                    heappush((other_v, w - x_values[edge]))
                heappush(other_v, (other_v, x_values[e]))

        # A and the remaining node v define a "cut of the phase". There is a
        # minimum cut of the original graph that is also a cut of the phase.
        # Due to contractions in earlier phases, v may in fact represent
        # multiple nodes in the original graph.
        # smallest value is heap(0)
        vert_weight = h[0]
        v = vert_weight[0]
        w = vert_weight[1]
        w = -w
        if w < cut_value:
            cut_value = w
            best_phase = i
        # Contract v and the last node added to A.
        contractions.append((u, v))
        for w, e in G[v].items():
            if w != u:
                if w not in G[u]:
                    G.add_edge(u, w, weight=e["weight"])
                else:
                    G[u][w]["weight"] += e["weight"]
        G.remove_node(v)

    # Recover the optimal partitioning from the contractions.
    # G = nx.Graph(islice(contractions, best_phase))
    # v = contractions[best_phase][1]
    # G.add_node(v)
    # reachable = set(nx.single_source_shortest_path_length(G, v))
    # partition = (list(reachable), list(nodes - reachable))

    return cut_value, partition