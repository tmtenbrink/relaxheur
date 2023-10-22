import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import gurobipy as gp
from enum import Enum, auto
from uuid import uuid4

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


def extended_formulation(n: int, graph_l: np.ndarray):
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

    cut_out_r = cut_out(0, n)
    cut_in_r = cut_in(0, n)

    for v in range(0, n):
        cut_idx = get_cut_idx(v, n)
        cut_char: gp.MVar = char_vec[cut_idx]
        model.addConstr(cut_char.sum() == 2, name=f"x(delta({v}))==2")

        # skip r
        if v == 0:
            continue

        s = v
        f_s: gp.MVar = flow[s - 1]
        fs_r_out: gp.MVar = f_s[cut_out_r]
        fs_r_in: gp.MVar = f_s[cut_in_r]
        model.addConstr(fs_r_out.sum() - fs_r_in.sum() >= 2, name=f"f_{s}(delta_out(r))-f_{s}(delta_in(r))>=2")

        for other_v in range(0, n):
            # skip r and s
            if other_v == 0 or other_v == s:
                continue

            cut_out_v = cut_out(other_v, n)
            cut_in_v = cut_in(other_v, n)

            fs_v_out: gp.MVar = f_s[cut_out_v]
            fs_v_in: gp.MVar = f_s[cut_in_v]

            model.addConstr(fs_v_out.sum() - fs_v_in.sum() == 0,
                            name=f"f_{s}(delta_out({other_v}))-f_{s}(delta_in({other_v}))==0")

        for m in range(m_edges):
            x_edge: gp.MVar = char_vec[m]
            f_s_arc_1: gp.MVar = f_s[m * 2]
            f_s_arc_2: gp.MVar = f_s[m * 2 + 1]

            model.addConstr(f_s_arc_1 - x_edge <= 0, name=f"f_{s}(arc {m * 2})<=x({m})")
            model.addConstr(f_s_arc_2 - x_edge <= 0, name=f"f_{s}(arc {m * 2 + 1})<=x({m})")

    return TSPModel(model, n, Formulation.EXTENDED, char_vec)


def cutting_plane_model(n: int, graph_l: np.ndarray):
    model = gp.Model("cutting plane")
    model.setParam('LogToConsole', 0)
    m_edges = (n * (n - 1)) // 2

    char_vec: gp.MVar = model.addMVar(shape=(m_edges,), lb=0, name='char_')

    z = (char_vec * graph_l).sum()
    model.setObjective(z, gp.GRB.MINIMIZE)

    return TSPModel(model, n, Formulation.CUTTING_PLANE)


class Formulation(Enum):
    EXTENDED = auto()
    CUTTING_PLANE = auto()


class TSPModel:
    formulation: Formulation
    is_relaxed = False
    model: gp.Model
    char_vec: Optional[gp.MVar]
    n: int

    def __init__(self, model, n: int, formulation: Formulation, char_vec: Optional[gp.MVar] = None):
        self.formulation = formulation
        self.is_relaxed = False
        self.model = model
        self.n = n

        self.char_vec = char_vec

    def optimize(self):
        if self.formulation == Formulation.EXTENDED:
            self.model.optimize()
            obj_value = self.model.ObjVal
            if self.is_relaxed:
                print(f"extended relaxation: {obj_value}")
            else:
                print(f"integer optimal: {int(obj_value)}")
        else:
            obj_value = optimize_cut_model(self)
            print(f"cutting plane relaxation: {obj_value}")

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
        if self.char_vec is None:
            var_list = []
            x_values = np.zeros(m_edges)
            for e in range(m_edges):
                char_e = self.model.getVarByName(f"char_[{e}]")
                var_list.append(char_e)
                x_values[e] = char_e.X
            return gp.MVar.fromlist(var_list), x_values
        else:
            x_values = np.zeros(m_edges)
            for e, char_var in enumerate(self.char_vec.tolist()):
                x_values[e] = char_var.X

            return self.char_vec, x_values

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


def compute_cut_weight(cut: np.ndarray, x_values: np.ndarray):
    # TODO return the weight = x(delta(cut))
    return 0


def compute_min_cut(x_values: np.ndarray):
    raise NotImplementedError()

    # TODO implement min cut

    weight = compute_cut_weight(np.ndarray([0]), x_values)
    # TODO return delta(U) and the cut_weight
    return [0], weight


def separation(n: int, char_vec: gp.MVar, x_values: np.ndarray, model: gp.Model) -> bool:
    """Tests whether x is in P_subtour, if not it adds the inequality to the model"""
    # based on bounds we put in model we assume x >= 0
    epsilon = 0.0001

    num_cstr = len(model.getConstrs())

    for v in range(n):
        edges = get_cut_idx(v, n)
        if np.abs(np.sum(x_values[edges]) - 2) > epsilon:
            edge_vars = char_vec[edges]
            model.addConstr(edge_vars.sum() == 2, name=f"x(delta({v}))==2")
            return False

    delta_cut, cut_weight = compute_min_cut(x_values)
    epsilon = 0.0001
    if 2 - cut_weight > epsilon:
        cut_vars = char_vec[delta_cut]
        model.addConstr(cut_vars.sum() >= 2, name=f"x(delta(cut_{uuid4()}))>=2")
        return False

    return True


def optimize_cut_model(m: TSPModel):
    max_i = 100000
    i = 0
    invalid = True
    while invalid:
        if i > max_i:
            raise ValueError("Model could not be solved in time!")

        m.model.optimize()
        char_vec, x_values = m.char_vec_values()
        # this modifies the underlying model and adds constraints
        in_subtour = separation(m.n, char_vec, x_values, m.model)

        invalid = not in_subtour
        i += 1

    return m.model.ObjVal


def run():
    # inst_path = get_inst_path()
    inst_path = Path('tsp/burma14.dat')

    n, graph_l = parse_instance(inst_path)

    cut_model = cutting_plane_model(n, graph_l)
    cut_model.optimize()
    cut_model.print_sol()

    # model = extended_formulation(n, graph_l)
    #
    # m = ExtendedModel(model, n)
    # m.optimize()
    # m.print_sol()

    # m_relax = m.relax()
    # m_relax.optimize()
    # m_relax.print_sol()

    # for i, x in enumerate(relax_model.getVars()):
    #     print(x)

    # we fix r is the first vertex


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
