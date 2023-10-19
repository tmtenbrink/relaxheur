import argparse
from pathlib import Path
import numpy as np
import gurobipy as gp


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

    return n, graph


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


def extended_formulation(n: int, graph: np.ndarray):
    model = gp.Model()
    m_edges = (n * (n - 1)) // 2

    lengths = np.zeros(m_edges, dtype=int)
    row_num = 0
    for row in graph:
        sel_from_row = np.arange(row_num + 1, n, dtype=int)
        selected = row[sel_from_row]
        target_start = row_num * (2 * n - row_num - 1) // 2
        target_end = target_start + selected.size
        lengths[target_start:target_end] = selected
        row_num += 1

    # for i<j (i, j) is the (n-1 + n-2 + ... + n - i) + (j - i -1)th element
    # i.e. index is i(2n-i-1)/2 + (j-i-1)
    char_vec: gp.MVar = model.addMVar(shape=(m_edges,), vtype=gp.GRB.BINARY)
    num_vert_non_r = n - 1
    num_arcs = 2 * m_edges
    # rows are vertices without r
    # columns are the arcs, with same indexing as above
    # but where (i, j) and (j, i) follow each other (so edge_index * 2 for (i, j)
    # with still i<j
    flow = model.addMVar(shape=(num_vert_non_r, num_arcs), lb=0)
    z = (char_vec * lengths).sum()
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
        f_s: gp.MVar = flow[s-1]
        fs_r_out: gp.MVar = f_s[cut_out_r]
        fs_r_in: gp.MVar = f_s[cut_in_r]
        model.addConstr(fs_r_out.sum() - fs_r_in.sum() == 0, name=f"f_{s}(delta_out(r))-f_{s}(delta_in(r))==0")

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

            model.addConstr(x_edge - f_s_arc_1 <= 0, f"f_{s}(arc {m * 2})<=x({m})")
            model.addConstr(x_edge - f_s_arc_2 <= 0, f"f_{s}(arc {m * 2 + 1})<=x({m})")

    return model, char_vec


def run():
    # inst_path = get_inst_path()
    path = Path('tsp/simple.dat')

    n, graph = parse_instance(path)

    model, char_vec = extended_formulation(n, graph)

    model.optimize()

    for x in char_vec.tolist():
        print(x.X)

    # we fix r is the first vertex


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
