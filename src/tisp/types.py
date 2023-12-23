from dataclasses import dataclass
from enum import Enum, auto
from functools import total_ordering

import mip

# EDGE INDEX CONVENTION:
# The lower node index precedes the higher node index. Start with node 0. Edge 0 is then (0, 1). Edge 1 is (0, 2), ... 
# up to Edge n-2 is (0, n-1). Then Edge n-1 is (1, 1) [we already had (1, 0)!], Edge n is (1, 2), ...
# This gives m=n*(n-1)/2 edges.
# If the lower index node is i and the higher index node is j, the edge index is:
# (i * (2n - i - 1) / 2) + (j - i - 1)

# Adjacency matrix for all costs. costs[nd_a][nd_b] will give the cost of the edge (nd_a, nd_b). We assume the problem 
# to be symmetric, so the cost matrix is also symmetric.
Costs = list[list[float]]
# Costs indexed by edge, using our edge index convention.
EdgeCosts = list[float]

# (lower index, higher index)
Edge = tuple[int, int]

# vertex_edges[v] gives the list of edges that contain the vertex v
VertexEdges = list[list[int]]
EdgesByIndex = dict[int, Edge]

# values of each edge (1 if included, 0 if included), using edge index convention
EdgeValues = list[float]

BestNeighbors = list[list[int]]

######### Tour

# begins with zero, followed by lowest index neighbor
Tour = list[int]


########## LP

class Formulation(Enum):
    EXTENDED = auto()
    CUTTING_PLANE = auto()

# n, m_edges, vertex_edges, edge_costs
LPConstants = tuple[int, int, VertexEdges, EdgeCosts, EdgesByIndex]

LPModel = tuple[mip.Model, Formulation, LPConstants]

EdgeValue = tuple[float, list[tuple[int, int]]]

######### Heuristic


HeurConstants = tuple[int, Costs, BestNeighbors]


######### Branch & Bound

PseudoList = tuple[list[float], list[float], list[int]]


HeurCosts = list[list[float]]
SubproblemState = tuple[LPModel, list[Edge], list[Edge], int, bool, float, HeurCosts]


Timer = tuple[int, int]

    
@total_ordering
@dataclass
class Subproblem:
    parent_lb: float
    state: SubproblemState

    # we define lt method so that they can be compared by heapsort
    def __lt__(self, other):
        return self.parent_lb < other.parent_lb
    


