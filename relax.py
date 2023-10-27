import argparse
from pathlib import Path
from typing import Optional
from enum import Enum, auto
from time import perf_counter_ns

import numpy as np
import gurobipy as gp


def parse_line(ln: str):
    return list(map(lambda i: int(i), ln.rstrip().split(" ")))


def get_inst_path():
    parser = argparse.ArgumentParser(description="Process some integers.")

    parser.add_argument("inst_name", help="Filename of instance")

    args = parser.parse_args()

    return Path(args.inst_name)


def parse_instance(inst_path: Path):
    with open(inst_path, "r") as f:
        lines = f.readlines()

    line_0 = parse_line(lines[0])
    n = line_0[0]
    all_lines = list(map(lambda ln: parse_line(ln), lines[1:]))
    graph = np.array(all_lines)

    m_edges = (n * (n - 1)) // 2

    graph_l = np.zeros(m_edges, dtype=int)
    row_num = 0
    for row in graph:
        sel_from_row = np.arange(row_num + 1, n, dtype=int)
        selected = row[sel_from_row]
        target_start = row_num * (2 * n - row_num - 1) // 2
        target_end = target_start + selected.size
        graph_l[target_start:target_end] = selected
        row_num += 1

    return n, graph_l


def get_edge_names(n: int):
    names = {}
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            names[index] = f"({i}, {j})"

            index += 1

    return names


def edge_idx(lower_i, higher_j, n: int):
    """lower_i and higher_j must be of same dimension or one of them must be a scalar."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)


def get_vert_indexes_for_cut(vertices: np.ndarray, n: int):
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

    return vertices_repeated, indexes_repeated, lt_mask, gt_mask


def get_cut_idx_arr(vertices: np.ndarray, n: int):
    """Vertices should be 1D array."""
    vertices_repeated, indexes_repeated, lt_mask, gt_mask = get_vert_indexes_for_cut(
        vertices, n
    )

    # we can now use edge_idx as before
    vertices_repeated[gt_mask] = edge_idx(
        vertices_repeated[gt_mask], indexes_repeated[gt_mask], n
    )
    vertices_repeated[lt_mask] = edge_idx(
        indexes_repeated[lt_mask], vertices_repeated[lt_mask], n
    )

    return vertices_repeated


def cut_out_arr(vertices: np.ndarray, n: int):
    vertices_repeated, indexes_repeated, lt_mask, gt_mask = get_vert_indexes_for_cut(
        vertices, n
    )

    # we can now use edge_idx as before
    vertices_repeated[gt_mask] = (
        edge_idx(vertices_repeated[gt_mask], indexes_repeated[gt_mask], n) * 2 + 1
    )
    vertices_repeated[lt_mask] = (
        edge_idx(indexes_repeated[lt_mask], vertices_repeated[lt_mask], n) * 2
    )

    return vertices_repeated


def cut_in_arr(vertices: np.ndarray, n: int):
    vertices_repeated, indexes_repeated, lt_mask, gt_mask = get_vert_indexes_for_cut(
        vertices, n
    )

    # we can now use edge_idx as before
    vertices_repeated[gt_mask] = (
        edge_idx(vertices_repeated[gt_mask], indexes_repeated[gt_mask], n) * 2
    )
    vertices_repeated[lt_mask] = (
        edge_idx(indexes_repeated[lt_mask], vertices_repeated[lt_mask], n) * 2 + 1
    )

    return vertices_repeated


def extended_formulation(n: int, graph_l: np.ndarray):
    print("Building extended formulation...")
    start_build = perf_counter_ns()
    model = gp.Model("extended")
    model.setParam("LogToConsole", 0)
    m_edges = (n * (n - 1)) // 2

    # for i<j (i, j) is the (n-1 + n-2 + ... + n - i) + (j - i -1)th element
    # i.e. index is i(2n-i-1)/2 + (j-i-1)
    char_vec: gp.MVar = model.addMVar(
        shape=(m_edges,), vtype=gp.GRB.BINARY, lb=0, name="char_"
    )
    num_vert_non_r = n - 1
    num_arcs = 2 * m_edges
    # rows are vertices without r
    # columns are the arcs, with same indexing as above
    # but where (i, j) and (j, i) follow each other (so edge_index * 2 for (i, j)
    # with still i<j
    flow = model.addMVar(shape=(num_vert_non_r, num_arcs), lb=0)
    z = (char_vec * graph_l).sum()
    model.setObjective(z, gp.GRB.MINIMIZE)

    # each column is all the vertices other than the vertex corresponding to column
    # so (n-1), n array
    cut_arr = get_cut_idx_arr(np.arange(n), n)
    repeated_vars = char_vec[cut_arr].reshape((n - 1, n))
    model.addConstr(repeated_vars.sum(axis=0) == 2, name=f"x(delta(v))==2")
    # we take transpose so we can easily access the rows
    cuts_in = cut_in_arr(np.arange(0, n), n).T
    cuts_out = cut_out_arr(np.arange(0, n), n).T

    cut_in_r = cuts_in[0]
    cut_out_r = cuts_out[0]
    fs_r_out: gp.MVar = flow[:, cut_out_r]
    fs_r_in: gp.MVar = flow[:, cut_in_r]
    model.addConstr(
        fs_r_out.sum(axis=1) - fs_r_in.sum(axis=1) >= 2,
        name=f"f_s(delta_out(r))-f_s(delta_in(r))>=2",
    )

    for s in range(1, n):
        f_s: gp.MVar = flow[s - 1]
        # we make sure we get arrays correspond to the right index of the char. vector
        arcs_one_direction = np.arange(m_edges) * 2
        arcs_other_direction = np.arange(m_edges) * 2 + 1
        model.addConstr(f_s[arcs_one_direction] - char_vec <= 0, name=f"f_1")
        model.addConstr(f_s[arcs_other_direction] - char_vec <= 0, name=f"f_2")

        for other_v in range(0, n):
            # skip r and s
            if other_v == 0 or other_v == s:
                continue

            cut_out_v = cuts_out[other_v]
            cut_in_v = cuts_in[other_v]

            fs_v_out: gp.MVar = f_s[cut_out_v]
            fs_v_in: gp.MVar = f_s[cut_in_v]

            model.addConstr(
                fs_v_out.sum() - fs_v_in.sum() == 0,
                name=f"f_{s}(delta_out({other_v}))-f_{s}(delta_in({other_v}))==0",
            )
    # We run update to make sure the timer works
    model.update()
    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return TSPModel(model, n, Formulation.EXTENDED)


def cutting_plane_model(n: int, graph_l: np.ndarray):
    print("Setting up cutting plane model...")
    start_build = perf_counter_ns()
    model = gp.Model("cutting plane")
    # set to dual simplex
    model.setParam("Method", 1)
    model.setParam("LogToConsole", 0)
    m_edges = (n * (n - 1)) // 2

    char_vec: gp.MVar = model.addMVar(shape=(m_edges,), lb=0, name="char_")

    z = (char_vec * graph_l).sum()
    model.setObjective(z, gp.GRB.MINIMIZE)
    model.update()
    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return TSPModel(model, n, Formulation.CUTTING_PLANE, True)


class Formulation(Enum):
    EXTENDED = auto()
    CUTTING_PLANE = auto()


class TSPModel:
    formulation: Formulation
    is_relaxed = False
    model: gp.Model
    char_vec: Optional[gp.MVar]
    n: int

    def __init__(self, model, n: int, formulation: Formulation, is_relaxed=False):
        self.formulation = formulation
        self.is_relaxed = is_relaxed
        self.model = model
        self.n = n

    def optimize(self):
        start_opt = perf_counter_ns()
        if self.formulation == Formulation.EXTENDED:
            relaxed_text = "relaxed " if self.is_relaxed else ""
            print(f"Optimizing {relaxed_text}extended formulation...")
            self.model.optimize()
            obj_value = self.model.ObjVal
            if self.is_relaxed:
                print(f"\t- extended relaxation: {obj_value}")
            else:
                print(f"\t- integer optimal: {int(obj_value)}")
        else:
            print("Optimizing cutting plane relaxation...")
            obj_value = optimize_cut_model(self)
            print(f"\t- cutting plane relaxation: {obj_value}")
        opt_time = (perf_counter_ns() - start_opt) / 1e9
        print(f"Optimizing took {opt_time} s.")

    def relax(self):
        self.model.update()
        relax_model: gp.Model = self.model.copy().relax()
        self.model = relax_model
        self.is_relaxed = True
        return self

    def char_vec_values(self) -> tuple[gp.MVar, np.ndarray]:
        m_edges = (self.n * (self.n - 1)) // 2
        var_list = []
        x_values = np.zeros(m_edges)
        for e in range(m_edges):
            char_e = self.model.getVarByName(f"char_[{e}]")
            var_list.append(char_e)
            x_values[e] = char_e.X
        return gp.MVar.fromlist(var_list), x_values

    def print_sol(self):
        char_vec, x_values = self.char_vec_values()
        edge_names = get_edge_names(self.n)

        epsilon = 0.0000001
        print("Path:")
        for e_i in range(x_values.size):
            x_val = x_values[e_i]
            if np.abs(x_val) > epsilon:
                add_value = f" with value {x_val}" if self.is_relaxed else ""
                print(f"\tedge {edge_names[e_i]} in solution{add_value}")

        return char_vec, x_values


def adjacency_matrix_from_vector(x_values: np.ndarray, n: int):
    """Create an adjacency matrix from x_values"""
    m_edges = (n * (n - 1)) // 2
    if len(x_values) != (n * (n - 1)) // 2:
        raise ValueError("x_values does not have the correct length")

    # adjacency_matrix = [[[0, [(i, j)]] for i in range(n)] for j in range(n)]
    adjacency_matrix = np.zeros((n, n))
    matr_edges = np.zeros((n, n, m_edges), dtype=int)

    for i in range(n):
        for j in range(i + 1, n):
            edge = edge_idx(i, j, n)
            matr_edges[i][j][edge] = 1
            matr_edges[j][i][edge] = 1

    return adjacency_matrix, matr_edges


def sw_minimum_cut_phase(graph_weights, graph_edges, a):
    graph_n = graph_weights.shape[0]
    A = [a]
    cut_edges = np.array([])
    max_cut_weight = -1

    while len(A) < graph_n:
        max_cut_weight = -1
        u = 0

        for v in range(graph_n):
            if v not in A:
                cut_weight = np.sum(graph_weights[v, A])
                if cut_weight > max_cut_weight:
                    max_cut_weight = cut_weight
                    u = v
                    cut_edges = np.sum(graph_edges[v, A], axis=0)

        A.append(u)

    s = min(A[-1], A[-2])
    t = max(A[-1], A[-2])

    return s, t, max_cut_weight, cut_edges


def sw_minimum_cut(graph_weights, graph_edges):
    """Find the minimum cut of a graph using the Stoer-Wagner algorithm."""
    n = graph_weights.shape[0]

    min_cut = float("inf")
    contractions = []
    phase = 0
    min_cut_edges = np.array([])

    while n > 1:
        a = 0  # Any vertex from V
        s, t, cut_weight, cut_edges = sw_minimum_cut_phase(graph_weights[:n, :n], graph_edges[:n, :n, :], a)
        if cut_weight < min_cut:
            min_cut = cut_weight
            min_cut_edges = cut_edges

        # Merge vertices s and t
        contractions.append((s, t))
        for i in range(n):
            if i != t:
                graph_weights[s, i] = graph_weights[s, i] + graph_weights[t, i]
                graph_edges[s, i] = graph_edges[s, i] + graph_edges[t, i]

        for i in range(t, n - 1):
            for j in range(n):
                graph_weights[i, j] = graph_weights[i + 1, j]
                graph_weights[j, i] = graph_weights[i][j]
                graph_edges[i, j] = graph_edges[i + 1, j]
                graph_edges[j, i] = graph_edges[i + 1, j]

        n -= 1
        phase += 1

    return min_cut, min_cut_edges


def compute_min_cut(x_values: np.ndarray, n: int) -> tuple[int, np.ndarray]:
    """Compute min cut using Stoer–Wagner algorithm. x_values should be 1D array with our edge
    index convention.
    """
    # we compute an adjacency matrix to more easily perform the stoer-wagner algorithm
    graph, graph_edges = adjacency_matrix_from_vector(x_values, n)
    min_cut, min_cut_edges = sw_minimum_cut(graph, graph_edges)
    # we turn the list of edges back into array using our edge index convention
    cut_edge_idx = np.flatnonzero(min_cut_edges)

    return min_cut, cut_edge_idx


def separation(
    n: int, char_vec: gp.MVar, x_values: np.ndarray, model: gp.Model
) -> bool:
    """Tests whether x is in P_subtour, if not it adds the violated inequality to the model.
    1) The constraint >=0 is already defined in the char_vector.
    2) We check if x(delta(v))=2, for all v, with some tolerance epsilon.
    3) We check x(delta(U))>= 2, for all U, by finding the minimum cut and checking if it
    is larger than 2 (with tolerance epsilon). Note that if the minimimum cut is larger
    than 2, we know that all cuts are larger than 2.

    It is easy to see that constraints 1 and 2 are checked in polynomial time. Constraint 3
    has exponentially many inequalities (as there are exponentially many U), but can be checked
    in polynomial time since the min-cut can be found in polynomial time by the Stoer-Wagner algorithm.

    Therefore, our separation algorithm is polynomial time.
    """

    # based on bounds we put in model we assume x >= 0
    epsilon = 0.0000001

    cuts = get_cut_idx_arr(np.arange(n), n)
    for v in range(n):
        # the columns are the vertices
        edges = cuts[:, v]
        if np.abs(np.sum(x_values[edges]) - 2) > epsilon:
            edge_vars = char_vec[edges]
            model.addConstr(edge_vars.sum() == 2, name=f"x(delta({v}))==2")
            return False

    cut_weight, min_cut_edges = compute_min_cut(x_values, n)
    if cut_weight < 2 - epsilon:
        subtour_vars = char_vec[min_cut_edges]
        model.addConstr(subtour_vars.sum() >= 2, name=f"x(delta(cut)>=2")
        return False

    return True


def optimize_cut_model(m: TSPModel):
    max_i = 100000
    i = 0
    invalid = True

    while invalid:
        if i > max_i:
            print("Taking a long time...")

        m.model.optimize()
        char_vec, x_values = m.char_vec_values()
        # this modifies the underlying model and adds constraints
        in_subtour = separation(m.n, char_vec, x_values, m.model)

        invalid = not in_subtour
        i += 1

    return m.model.ObjVal


def run():
    # inst_path = get_inst_path()
    inst_path = Path('tsp/gr48.dat')

    n, graph_l = parse_instance(inst_path)

    # make empty model to startup Gurobi
    _ = gp.Model()

    cut_model = cutting_plane_model(n, graph_l)
    cut_model.optimize()
    # cut_model.print_sol()

    model = extended_formulation(n, graph_l)
    model.optimize()
    # model.print_sol()

    m_relax = model.relax()
    m_relax.optimize()
    # m_relax.print_sol()


if __name__ == "__main__":
    run()
