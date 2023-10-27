import numpy as np


def min_cut_phase(edge_weights: np.ndarray, w, a, edge_tuples_arr: np.ndarray, n: int, cut_arr: np.ndarray):
    A = {a}
    cut = cut_arr[a]


