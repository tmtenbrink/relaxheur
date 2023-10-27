import argparse
from pathlib import Path
from typing import Optional
from enum import Enum, auto
from uuid import uuid4
from time import perf_counter_ns

import numpy as np
import gurobipy as gp
import networkx as nx


def parse_line(ln: str):
    return list(map(lambda i: int(i), ln.rstrip().split(' ')))


def get_inst_path():
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('inst_name', help='Filename of instance')

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


def vertex_from_edge(edge: int):
    # i, j
    # I = i * (2*n - i - 1) / 2 + (j - i - 1)
    return


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


def get_vert_indexes_for_cut(vertices: np.ndarray, n: int):
    vertices_len = vertices.size
    # each column is the vertex repeated
    vertices_repeated = np.tile(vertices, (n - 1, 1))
    # each column is the index from 0 to n-1
    indexes_repeated = np.repeat(np.arange(n - 1), vertices_len).reshape((n - 1, vertices_len))
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
    vertices_repeated, indexes_repeated, lt_mask, gt_mask = get_vert_indexes_for_cut(vertices, n)

    # we can now use edge_idx as before
    vertices_repeated[gt_mask] = edge_idx(vertices_repeated[gt_mask], indexes_repeated[gt_mask], n)
    vertices_repeated[lt_mask] = edge_idx(indexes_repeated[lt_mask], vertices_repeated[lt_mask], n)

    return vertices_repeated


def edge_idx_to_arc_idx(edge_idx_arr: np.ndarray):
    """indexes must be a 1D array"""
    num_edges = edge_idx_arr.size
    arc_idx_arr = np.zeros(2 * num_edges, dtype=int)
    split_indexes = np.arange(start=0, stop=num_edges, dtype=int)
    base_idx = edge_idx_arr * 2
    reverse_idx = base_idx + 1
    arc_idx_arr[split_indexes] = base_idx
    arc_idx_arr[split_indexes + 1] = reverse_idx

    return arc_idx_arr


def cut_out(vertex_i: int, n: int):
    # (u, v) means u is tail
    # so cut_out is all (u, v) for v not equal to u
    less_idx, more_idx = get_lt_and_gt_cut(vertex_i, n)

    arc_cut_size = (n - 1)
    cut_out_arr = np.zeros(arc_cut_size, dtype=int)
    # for edges where other vertex has lower index, we need the second arc, because vertex_i must be first
    cut_out_arr[:vertex_i] = less_idx * 2 + 1
    # for edges where vertex_i has lower index, we need the first arc, because vertex_i must be first
    cut_out_arr[vertex_i:] = more_idx * 2

    return cut_out_arr


def cut_in(vertex_i: int, n: int):
    # (u, v) means v is head
    # so cut_in is all (u, v) for u not equal to v
    less_idx, more_idx = get_lt_and_gt_cut(vertex_i, n)

    arc_cut_size = (n - 1)
    cut_in_arr = np.zeros(arc_cut_size, dtype=int)
    # for edges where other vertex has lower index, we need the first, because vertex_i must be second
    cut_in_arr[:vertex_i] = less_idx * 2
    # for edges where vertex_i has lower index, we need the second arc, because vertex_i must be second
    cut_in_arr[vertex_i:] = more_idx * 2 + 1

    return cut_in_arr


def cut_out_arr(vertices: np.ndarray, n: int):
    vertices_repeated, indexes_repeated, lt_mask, gt_mask = get_vert_indexes_for_cut(vertices, n)

    # we can now use edge_idx as before
    vertices_repeated[gt_mask] = edge_idx(vertices_repeated[gt_mask], indexes_repeated[gt_mask], n) * 2 + 1
    vertices_repeated[lt_mask] = edge_idx(indexes_repeated[lt_mask], vertices_repeated[lt_mask], n) * 2

    return vertices_repeated


def cut_in_arr(vertices: np.ndarray, n: int):
    vertices_repeated, indexes_repeated, lt_mask, gt_mask = get_vert_indexes_for_cut(vertices, n)

    # we can now use edge_idx as before
    vertices_repeated[gt_mask] = edge_idx(vertices_repeated[gt_mask], indexes_repeated[gt_mask], n) * 2
    vertices_repeated[lt_mask] = edge_idx(indexes_repeated[lt_mask], vertices_repeated[lt_mask], n) * 2 + 1

    return vertices_repeated


def get_edge_tpls_arr(n: int):
    tpls = []
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            tpls.append((i, j, 0.0))

            index += 1

    return np.array(tpls)


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


def extended_formulation(n: int, graph_l: np.ndarray):
    print("Building extended formulation...")
    start_build = perf_counter_ns()
    model = gp.Model("extended")
    model.setParam('LogToConsole', 0)
    m_edges = (n * (n - 1)) // 2

    # for i<j (i, j) is the (n-1 + n-2 + ... + n - i) + (j - i -1)th element
    # i.e. index is i(2n-i-1)/2 + (j-i-1)
    char_vec: gp.MVar = model.addMVar(shape=(m_edges,), vtype=gp.GRB.BINARY, lb=0, name='char_')
    num_vert_non_r = n - 1
    num_arcs = 2 * m_edges
    # rows are vertices without r
    # columns are the arcs, with same indexing as above
    # but where (i, j) and (j, i) follow each other (so edge_index * 2 for (i, j)
    # with still i<j
    flow = model.addMVar(shape=(num_vert_non_r, num_arcs), lb=0)
    z = (char_vec * graph_l).sum()
    model.setObjective(z, gp.GRB.MINIMIZE)

    # cut_out_r = cut_out(0, n)
    # cut_in_r = cut_in(0, n)
    build_i = 0
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
    model.addConstr(fs_r_out.sum(axis=1) - fs_r_in.sum(axis=1) >= 2, name=f"f_s(delta_out(r))-f_s(delta_in(r))>=2")

    # for v in range(n):
    #     print(f"in i {cuts_in[v]}")
    #     print(f"in {cut_in(v, n)}")
    #
    #     print(f"out i {cuts_out[v]}")
    #     print(f"out {cut_in(v, n)}")

    for s in range(1, n):
        f_s: gp.MVar = flow[s - 1]
        for other_v in range(0, n):
            # skip r and s
            if other_v == 0 or other_v == s:
                continue

            cut_out_v = cuts_out[other_v]
            cut_in_v = cuts_in[other_v]
            # cut_out_v = cut_out(other_v, n)
            # cut_in_v = cut_in(other_v, n)

            fs_v_out: gp.MVar = f_s[cut_out_v]
            fs_v_in: gp.MVar = f_s[cut_in_v]

            model.addConstr(fs_v_out.sum() - fs_v_in.sum() == 0,
                            name=f"f_{s}(delta_out({other_v}))-f_{s}(delta_in({other_v}))==0")
        # we make sure we get arrays correspond to the right index of the char. vector
        arcs_one_direction = np.arange(m_edges) * 2
        arcs_other_direction = np.arange(m_edges) * 2 + 1
        model.addConstr(f_s[arcs_one_direction] - char_vec <= 0, name=f"f_1")
        model.addConstr(f_s[arcs_other_direction] - char_vec <= 0, name=f"f_2")
    model.update()
    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return TSPModel(model, n, Formulation.EXTENDED)


def cutting_plane_model(n: int, graph_l: np.ndarray):
    print("Setting up cutting plane model...")
    start_build = perf_counter_ns()
    model = gp.Model("cutting plane")
    model.setParam('LogToConsole', 0)
    m_edges = (n * (n - 1)) // 2

    char_vec: gp.MVar = model.addMVar(shape=(m_edges,), lb=0, name='char_')

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
                print(f"extended relaxation: {obj_value}")
            else:
                print(f"integer optimal: {int(obj_value)}")
        else:
            print("Optimizing cutting plane relaxation...")
            obj_value = optimize_cut_model(self)
            print(f"cutting plane relaxation: {obj_value}")
        opt_time = (perf_counter_ns() - start_opt) / 1e9
        print(f"Optimizing took {opt_time} s.")

    def relax(self):
        self.model.update()
        relax_model: gp.Model = self.model.copy().relax()
        # model_vars = relax_model.getVars()
        # for v in model_vars:
        #     print(v.VType)
        # model_cstrs = relax_model.getConstrs()
        # for c in model_cstrs:
        #     print(c.ConstrName)
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

        epsilon = 0.0001
        print("Path:")
        for e_i in range(x_values.size):
            x_val = x_values[e_i]
            if np.abs(x_val) > epsilon:
                add_value = f" with value {x_val}" if self.is_relaxed else ""
                print(f"\tedge {edge_names[e_i]} in solution{add_value}")

        return char_vec, x_values


# EDGE_24 = get_edge_names(24)
#
# def edge_name_str(x_values):
#     epsilon = 0.0001
#     print("Path:")
#     for e_i in range(x_values.size):
#         x_val = x_values[e_i]
#         if x_val > epsilon:
#             add_value = f" with value {x_val}"
#             print(f"\tedge {EDGE_24[e_i]} in solution{add_value}")

def compute_min_cut(x_values: np.ndarray, n: int, base_graph: nx.Graph, edge_tuples_arr: np.ndarray) -> tuple[
    int, np.ndarray]:
    """Compute min cut using Stoer–Wagner algorithm using Networkx. x_values should be 1D array with our edge
    index convention. edge_tuples_arr should be of the same length, but here each element is a 3-element array
    (so an m*3 array) that is (i, j, w) with i the lower vertex index, j the higher and w to be used for weight."""
    # We make a copy so we do not modify the base graph
    weighted_graph = base_graph.copy()
    # We copy the array with weights zero
    edge_tuples_arr = edge_tuples_arr.copy()
    edge_tuples_arr[:, 2] = x_values
    weighted_graph.add_weighted_edges_from(edge_tuples_arr)

    cut_value, partition = nx.stoer_wagner(weighted_graph)
    cut_edges = partition_to_edge_idx(n, partition)
    # edge_name = edge_name_str(x_values)
    return cut_value, cut_edges


def separation(n: int, char_vec: gp.MVar, x_values: np.ndarray, model: gp.Model, base_graph: nx.Graph,
               edge_tuples_arr: np.ndarray) -> bool:
    """Tests whether x is in P_subtour, if not it adds the inequality to the model"""
    # based on bounds we put in model we assume x >= 0
    epsilon = 0.0001

    cstrs = model.getConstrs()
    num_cstr = len(cstrs)
    # print(num_cstr)

    for v in range(n):
        edges = get_cut_idx(v, n)
        if np.abs(np.sum(x_values[edges]) - 2) > epsilon:
            edge_vars = char_vec[edges]
            model.addConstr(edge_vars.sum() == 2, name=f"x(delta({v}))==2")
            return False

    cut_weight, min_cut_edges = compute_min_cut(x_values, n, base_graph, edge_tuples_arr)
    if cut_weight < 2:
        subtour_vars = char_vec[min_cut_edges]
        model.addConstr(subtour_vars.sum() >= 2, name=f"x(delta(cut_{uuid4()}))>=2")
        return False

    return True


def optimize_cut_model(m: TSPModel):
    max_i = 100000
    i = 0
    invalid = True

    edge_tpls_arr = get_edge_tpls_arr(m.n)
    base_graph = nx.complete_graph(m.n)
    while invalid:
        if i > max_i:
            raise ValueError("Model could not be solved in time!")

        m.model.optimize()
        char_vec, x_values = m.char_vec_values()
        # this modifies the underlying model and adds constraints
        in_subtour = separation(m.n, char_vec, x_values, m.model, base_graph, edge_tpls_arr)

        invalid = not in_subtour
        i += 1

    return m.model.ObjVal


def run():
    inst_path = get_inst_path()
    # inst_path = Path('tsp/bays29.dat')

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

    # we fix r is the first vertex


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
