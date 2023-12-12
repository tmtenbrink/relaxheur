import argparse
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from time import perf_counter_ns
from typing import Optional, Sequence, SupportsFloat
import gurobipy as gp

Costs = list[list[float]]
EdgeCosts = list[float]

def parse_line(ln: str) -> list[float]:
    return list(map(lambda i: float(i), ln.rstrip().split(" ")))


def get_inst_path():
    parser = argparse.ArgumentParser(description="Process some integers.")

    parser.add_argument("inst_name", help="Filename of instance")

    args = parser.parse_args()

    return Path(args.inst_name)


def parse_as_adj_matrix(inst_path: Path) -> tuple[int, Costs]:
    with open(inst_path, "r") as f:
        lines = f.readlines()

    line_0 = parse_line(lines[0])
    n = int(line_0[0])
    adj_matrix = list(map(lambda ln: parse_line(ln), lines[1:]))

    return n, adj_matrix


EdgeValue = tuple[float, list[tuple[int, int]]]


def adjacency_matrix_from_vector(x_values: list[float], n: int) -> list[list[EdgeValue]]:
    """Create an adjacency matrix from x_values"""
    if len(x_values) != (n * (n - 1)) // 2:
        raise ValueError("x_values does not have the correct length")

    adjacency_matrix = [[(0.0, [(i, j)]) for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            idx = edge_idx(i, j, n)
            adjacency_matrix[i][j] = (x_values[idx], [(i, j)])
            adjacency_matrix[j][i] = adjacency_matrix[i][j]

    return adjacency_matrix


def compute_edge_costs(costs: Costs) -> EdgeCosts:
    n = len(costs)
    m_edges = (n * (n - 1)) // 2
    edge_costs = [0.0]*m_edges

    for row_num, row in enumerate(costs):
        target_start = row_num * (2 * n - row_num - 1) // 2
        target_end = target_start + (n - row_num  - 1)
        for i, e in enumerate(range(target_start, target_end)):
            edge_costs[e] = row[row_num+1+i]

    return edge_costs



def get_edge_names(n: int):
    names = {}
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            names[index] = f"({i}, {j})"

            index += 1

    return names


def edge_idx(lower_i: int, higher_j: int, n: int):
    """Returns the index of the edge using our edge-index convention."""
    return lower_i * (2 * n - lower_i - 1) // 2 + (higher_j - lower_i - 1)
    


def get_edge_idxs_for_v(vertex: int, n: int):
    """Return the indexes of the """
    lower_edges = [edge_idx(other_v, vertex, n) for other_v in range(0, vertex)]
    higher_edges = [edge_idx(vertex, other_v, n) for other_v in range(vertex+1, n)]

    return lower_edges + higher_edges


def get_arc_idxs_for_v(vertex: int, n: int, incoming: bool):
    lower_add = 0 if incoming else 1
    higher_add = 1 - lower_add

    lower_edges = [2*edge_idx(other_v, vertex, n)+lower_add for other_v in range(0, vertex)]
    higher_edges = [2*edge_idx(vertex, other_v, n)+higher_add for other_v in range(vertex+1, n)]

    return lower_edges + higher_edges

def get_in_arc_idxs_for_v(vertex: int, n: int):
    """Return the indexes of the incoming arcs for a vertex, i.e. where v is the head (so second vertex in the ordered pair)."""
    return get_arc_idxs_for_v(vertex, n, True)


def get_out_arc_idxs_for_v(vertex: int, n: int):
    """Return the indexes of the incoming arcs for a vertex, i.e. where v is the head (so second vertex in the ordered pair)."""
    return get_arc_idxs_for_v(vertex, n, False)


def edge_idxs_for_all_v(n: int) -> list[list[int]]:
    return [get_edge_idxs_for_v(v, n) for v in range(n)]


class Formulation(Enum):
    EXTENDED = auto()
    CUTTING_PLANE = auto()


@dataclass
class Problem:
    n: int
    # number of edges
    m_edges: int
    # all edges containing a certain vertex (edges according to our index convention)
    all_cuts: list[list[int]]
    # list of coefficients for char_vec (i.e. [1]*m_edges)
    coeff: list[int]
    edge_costs: list[float]

    def all(self) -> tuple[int, int, list[list[int]], list[int]]:
        return self.n, self.m_edges, self.all_cuts, self.coeff



class UnknownVariableError(Exception):
    pass


class TSPModel:
    formulation: Formulation
    model: gp.Model
    char_vec: Optional[list[gp.Var]]
    p: Problem

    def __init__(self, model, n: int, edge_costs: EdgeCosts, formulation: Formulation):
        self.formulation = formulation
        self.model = model
        m_edges = m_edges = (n * (n - 1)) // 2
        all_cuts = edge_idxs_for_all_v(n)
        self.p = Problem(n, m_edges, all_cuts, [1]*len(all_cuts[0]), edge_costs)

    def optimize(self):
        start_opt = perf_counter_ns()
        if self.formulation == Formulation.EXTENDED:
            print(f"Optimizing integer extended formulation...")
            self.model.optimize()
            obj_value = self.model.ObjVal
            print(f"\t- integer optimal: {round(obj_value)}")
            relaxed_model = self.model.relax()
            relaxed_model.optimize()
            relaxed_obj_value = relaxed_model.ObjVal
            print(f"\t- extended relaxation: {relaxed_obj_value}")
        else:
            print("Optimizing cutting plane relaxation...")
            self.model = optimize_cut_model(self)
            obj_value = self.model.ObjVal
            print(f"\t- cutting plane relaxation: {obj_value}")
        opt_time = (perf_counter_ns() - start_opt) / 1e9
        print(f"Optimizing took {opt_time} s.")
        # self.print_sol()

    def relax(self):
        self.model.update()
        relax_model: gp.Model = self.model.copy().relax()
        self.model = relax_model
        self.is_relaxed = True
        return self

    def char_vec_values(self) -> tuple[list[gp.Var], list[float]]:
        m_edges = self.p.m_edges
        var_list = []
        x_values = [0.0]*m_edges
        for e in range(m_edges):
            char_e = self.model.getVarByName(f"char_{e}")
            if char_e is None:
                raise UnknownVariableError(f"Could not find variable char_{e}!")
            var_list.append(char_e)
            x_values[e] = char_e.X
        return var_list, x_values

    def print_sol(self):
        char_vec, x_values = self.char_vec_values()
        edge_names = get_edge_names(self.p.n)
        edge_costs = self.p.edge_costs

        epsilon = 0.0000001
        print("Path:")
        cost = 0
        for e_i in range(len(x_values)):
            x_val = x_values[e_i]
            if abs(x_val) > epsilon:
                add_value = f" with value {x_val}"
                print(f"\tedge {edge_names[e_i]} in solution{add_value}")
                cost += edge_costs[e_i]*x_val

        print(f"Cost: {cost}")
        return char_vec, x_values
    


def separation(
    p: Problem, char_vec: list[gp.Var], x_values: list[float], model: gp.Model
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

    n, _, cuts, coeff = p.all()
    
    for v in range(n):
        # the columns are the vertices
        edges = cuts[v]
        x_sum = 0
        for e in edges:
            x_sum += x_values[e]

        if abs(x_sum - 2) > epsilon:
            edge_vars = [char_vec[e] for e in edges]
            model.addConstr(gp.LinExpr(coeff, edge_vars) == 2, name=f"x(delta({v}))==2")
            return False

    cut_weight, min_cut_edges = compute_min_cut(x_values, n)
    if cut_weight < 2 - epsilon:
        subtour_vars = [char_vec[e] for e in min_cut_edges]
        coeff_subtour = [1]*len(subtour_vars)
        model.addConstr(gp.LinExpr(coeff_subtour, subtour_vars) >= 2, name=f"x(delta(cut)>=2")
        return False

    return True


def compute_min_cut(x_values: list[float], n: int) -> tuple[float, list[int]]:
    """Compute min cut using Stoer–Wagner algorithm. x_values should be 1D array with our edge
    index convention.
    """
    # we compute an adjacency matrix to more easily perform the stoer-wagner algorithm
    graph = adjacency_matrix_from_vector(x_values, n)
    min_cut, best_edge_list = sw_minimum_cut(graph)
    # we turn the list of edges back into array using our edge index convention
    cut_edges = [edge_idx(e[0], e[1], n) for e in best_edge_list]

    return min_cut, cut_edges



def sw_minimum_cut_phase(graph: list[list[EdgeValue]], a: int) -> tuple[int, int, float, list[tuple[int, int]]]:
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
        in_subtour = separation(m.p, char_vec, x_values, m.model)

        invalid = not in_subtour
        i += 1

    return m.model


def extended_formulation(n: int, edge_costs: EdgeCosts):
    model = gp.Model("extended")
    # don't log
    model.setParam("LogToConsole", 0)
    print("Building extended formulation...")
    start_build = perf_counter_ns()

    m_edges = (n * (n - 1)) // 2
    # for i<j (i, j) is the (n-1 + n-2 + ... + n - i) + (j - i -1)th element
    # i.e. index is i(2n-i-1)/2 + (j-i-1)
    char_vec = [model.addVar(lb=0, name=f"char_{e}", vtype=gp.GRB.BINARY) for e in range(m_edges)]
    num_arcs = 2 * m_edges
    # rows are vertices without r
    # columns are the arcs, with same indexing as above
    # but where (i, j) and (j, i) follow each other (so edge_index * 2 for (i, j)
    # with still i<j
    flow: list[list[gp.Var]] = [[model.addVar(lb=0, name=f"f_{s}({a})") for a in range(num_arcs)] for s in range(1, n)]
    z = gp.LinExpr(edge_costs, char_vec)
    model.setObjective(z, gp.GRB.MINIMIZE)

    # each column is all the vertices other than the vertex corresponding to column
    # so (n-1), n array
    cut_arr = edge_idxs_for_all_v(n)
    coeff_1 = [1]*len(cut_arr[0])
    for v in range(n):
        v_cut = cut_arr[v]
        cut_x = [char_vec[e] for e in v_cut]
        model.addConstr(gp.LinExpr(coeff_1, cut_x) == 2, name=f"x(delta({v}))==2")
    # we take transpose so we can easily access the rows
    cuts_in = [get_in_arc_idxs_for_v(v, n) for v in range(n)]
    cuts_out = [get_out_arc_idxs_for_v(v, n) for v in range(n)]

    cut_in_r = cuts_in[0]
    cut_out_r = cuts_out[0]
    coeff_r_1_in = [1]*len(cut_in_r)
    coeff_r_1_out = [1]*len(cut_out_r)
    
    for s in range(1, n):
        f_s = flow[s - 1]
        fs_r_out = [f_s[a] for a in cut_out_r]
        fs_r_in = [f_s[a] for a in cut_in_r]
        out_sum = gp.LinExpr(coeff_r_1_out, fs_r_out)
        in_sum = gp.LinExpr(coeff_r_1_in, fs_r_in)
        model.addConstr(
            out_sum - in_sum >= 2,
            name=f"f_{s}(delta_out(r))-f_{s}(delta_in(r))>=2",
        )


    for s in range(1, n):
        f_s = flow[s - 1]
        # we make sure we get arrays correspond to the right index of the char. vector
        for e in range(m_edges):
            a1 = 2*e
            a2 = 2*e + 1
            model.addConstr(f_s[a1] - char_vec[e] <= 0, name=f"f_{a1}<=x({a1})")
            model.addConstr(f_s[a2] - char_vec[e] <= 0, name=f"f_{a2}<=x({a2})")

        for other_v in range(0, n):
            # skip r and s
            if other_v == 0 or other_v == s:
                continue

            cut_out_v = cuts_out[other_v]
            cut_in_v = cuts_in[other_v]
            coeff_v_1_in = [1]*len(cut_in_v)
            coeff_v_1_out = [1]*len(cut_out_v)
            fs_v_out = [f_s[v] for v in cut_out_v]
            fs_v_in = [f_s[v] for v in cut_in_v]

            out_sum = gp.LinExpr(coeff_v_1_out, fs_v_out)
            in_sum = gp.LinExpr(coeff_v_1_in, fs_v_in)

            model.addConstr(
                out_sum - in_sum == 0,
                name=f"f_{s}(delta_out({other_v}))-f_{s}(delta_in({other_v}))==0",
            )
    # We run update to make sure the timer works
    model.update()
    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return TSPModel(model, n, edge_costs, Formulation.EXTENDED)



def cutting_plane_model(n: int, edge_costs: EdgeCosts):
    model = gp.Model("cutting plane")
    # set to dual simplex
    model.setParam("Method", 1)
    # don't log
    model.setParam("LogToConsole", 0)
    print("Setting up cutting plane model...")
    start_build = perf_counter_ns()
    
    m_edges = (n * (n - 1)) // 2
    char_vec = [model.addVar(lb=0, name=f"char_{e}") for e in range(m_edges)]
    z = gp.LinExpr(edge_costs, char_vec)
    model.setObjective(z, gp.GRB.MINIMIZE)
    
    model.update()
    build_time = (perf_counter_ns() - start_build) / 1e9
    print(f"Took {build_time} s.")
    return TSPModel(model, n, edge_costs, Formulation.CUTTING_PLANE)


def run():
    # inst_path = get_inst_path()
    inst_path = Path('tsp/gr24.dat')

    n, graph_l = parse_as_adj_matrix(inst_path)
    edge_costs = compute_edge_costs(graph_l)

    cut_model = cutting_plane_model(n, edge_costs)

    ext_model = extended_formulation(n, edge_costs)
    
    cut_model.optimize()
    ext_model.optimize()


if __name__ == "__main__":
    run()
