"""Very simple branch-and-bound for the knapsack problem.

PROBLEM:

    Given: a finite set N of items, a positive weight function 
           w: N -> R, a cost function c: N -> R, and a capacity U;

    Find : a subset S of N with w(S) <= U with maximum cost c(S).

SUBPROBLEMS:

    Subproblems are given by two disjoint subsets D and F of N. The
    set D contains all items that have to be in the solution and the
    set F contains all items that cannot be in the solution (i.e., D
    contains the variables that are fixed to 1 and F the variables
    that are fixed to 0). Items not in D nor F may or may not be in
    the solution.

UPPER AND LOWER BOUNDS FOR A SUBPROBLEM (D, F):

    Let (N, w, c, U) be an instance of the knapsack problem. Assume
    without loss of generality that N = { 1, ..., n } and that 

        c[1] / w[1] >= c[2] / w[2] >= ... >= c[n] / w[n].

    Write

        k = max { i : w[1] + ... + w[i] <= U }.

    A feasible solution for (N, w, c, U) is S = { 1, ..., k }, which
    gives us a lower bound of c(S). To get a lower bound for a
    subproblem (D, F), we apply this idea to the instance

    (*) (N - D - F, w, c, U - w(D)).

    If k = n, then w(N) <= U and c(N) is the optimal value of the
    instance (N, w, c, U). If not, then write

        x = (U - (w[1] + ... + w[k])) / w[k + 1].

    Theorem. The optimal value of the instance is at most

        c[1] + ... + c[k] + x * c[k + 1].

    (In other words, one orders the items by cost per unit of
    weight. Then one greedily selects all items starting from 1 while
    they fit in the knapsack. When the next item does not fit, one
    gets the maximum fraction of it that fits.)

    EXERCISE: Prove the theorem. (See Dantzig, "Discrete variables
    extremum problems", Operations Research 5 (1957) 266-277.)

    To get an upper bound for a subproblem (D, F), again we apply this
    idea to the instance (*).

AUTHOR:

    - Fernando M. de Oliveira Filho <fmario@gmail.com>

CREATED: <Thu Aug 30 13:40:30 2018 02:00>

"""

from random import randrange


class KnapsackInstance:
    def __init__(self, w, c, U):
        n = self.n = len(w)

        # We sort the weights and costs according to cost / weight
        # ratio.
        N = self.N = list(sorted(range(n), key = lambda i: -1.0 * c[i] / w[i]))
        self.w = [ w[i] for i in N ]
        self.c = [ c[i] for i in N ]
        self.U = U

        
def generate_random_instance(nitems):
    w = [ randrange(1, 101) for i in xrange(nitems) ]
    c = [ randrange(1, 5001) for i in xrange(nitems) ]
    return KnapsackInstance(w, c, 0.75 * sum(w))
        
    
class Subproblem:
    def __init__(self, fixed_one, fixed_zero, ub):
        self.fixed_one = set(fixed_one)
        self.fixed_zero = set(fixed_zero)
        self.ub = ub


class InfeasibleSubproblem(Exception):
    pass


def compute_bounds(inst, P):
    lb = 0
    total_weight = 0

    # First we add up the fixed items.
    for i in P.fixed_one:
        lb += inst.c[i]
        total_weight += inst.w[i]

    if total_weight > inst.U:
        raise InfeasibleSubproblem()

    # Then greedily fill up the rest to get a lower bound.
    items = []
    for i in xrange(inst.n):
        if i not in P.fixed_one and i not in P.fixed_zero:
            if total_weight + inst.w[i] <= inst.U:
                lb += inst.c[i]
                total_weight += inst.w[i]
                items.append(i)
            else:
                # First item not to fit.
                first_fail = i
                break
    else:
        # All items fit, so we have the optimal solution (upper bound
        # is equal to lower bound).
        return lb, lb, items

    # Compute the upper bound.
    ub = (lb + 1.0 * inst.c[first_fail]
          * (inst.U - total_weight) / inst.w[first_fail])

    return lb, ub, items


def do_branch_and_bound(inst):
    # Global upper and lower bounds, and best solution found so far.
    global_ub = sum(inst.c) + 1
    global_lb = -1
    best_solution = []

    # After processing a node, the new global upper bound is the
    # maximum of the upper bounds of active nodes and integer
    # nodes. Since integer nodes get out of the list of active nodes,
    # we keep the maximum upper bound of integer nodes in the
    # following variable.
    integer_node_bound = -1
    
    # Initialization.
    active_nodes = [ Subproblem([], [], global_ub) ]

    # Main loop.
    while active_nodes:
        # Select an active node to process.
        P = active_nodes[0]
        active_nodes = active_nodes[1:]

        # Process the node.
        try:
            lb, ub, items = compute_bounds(inst, P)
        except InfeasibleSubproblem:
            # Pruned by infeasibility.
            continue

        # Update global lower bound.
        if lb > global_lb:
            global_lb = lb
            best_solution = list(P.fixed_one) + items
            print 'Improved lower bound:', global_lb

        # Update global upper bound.
        if lb == ub and lb > integer_node_bound:
            integer_node_bound = lb
        
        if active_nodes:
            new_global_ub = max(ub, integer_node_bound,
                                max(P.ub for P in active_nodes))
        else:
            new_global_ub = max(ub, integer_node_bound)

        if new_global_ub < global_ub:
            global_ub = new_global_ub
            print 'Improved upper bound:', global_ub

        # Prune by bound?
        if ub < global_lb:
            continue
            
        # Prune by optimality?
        if lb == ub:
            continue

        # Select variable for split and perform the split.
        for i in xrange(inst.n):
            if i not in P.fixed_one and i not in P.fixed_zero:
                break
        else:
            raise RuntimeError('no variable to fix; this is a bug')

        Pl = Subproblem(list(P.fixed_one) + [ i ], P.fixed_zero, ub)
        Pr = Subproblem(P.fixed_one, list(P.fixed_zero) + [ i ], ub)
        active_nodes += [ Pl, Pr ]

    assert global_ub >= global_lb
        
    # Check that the solution is truly feasible.
    if sum(inst.w[i] for i in best_solution) > inst.U:
        raise RuntimeError('solution is infeasible; this is a bug')
        
    # Return optimal solution.
    return lb, best_solution


def run_gurobi(inst):
    """Runs Gurobi to solve the knapsack problem."""

    try:
        from gurobipy import Model, GRB
    except:
        raise RuntimeError('Gurobi not found')
    
    model = Model('')

    # Variables.
    x = model.addVars(range(inst.n), vtype = GRB.BINARY)

    # Add capacity constraint.
    lhs = 0
    for i in xrange(inst.n):
        lhs += inst.w[i] * x[i]

    model.addConstr(lhs <= inst.U)

    # Objective function.
    model.setObjective(sum(inst.c[i] * x[i] for i in xrange(inst.n)),
                       GRB.MAXIMIZE)

    model.optimize()
    

def main():
    inst = generate_random_instance(10)
    #w = [ 7, 3, 17, 35, 81 ]
    #c = [ 4416, 211, 849, 1077, 445 ]
    #inst = KnapsackInstance(w, c, 107)

    print 'Here is the output of the branch-and-bound method'
    opt, sol = do_branch_and_bound(inst)

    print '\nOptimal solution =', opt, ' items =', sol

    print '\nNow we run Gurobi...\n'

    try:
        run_gurobi(inst)
    except RuntimeError as re:
        print 'error:', re.message


main()
