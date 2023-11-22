import random

def dict_tour(tour: list[int]) -> dict[int, tuple[int, int]]:
    # adjacent edges in the tour
    d: dict[int, tuple[int, int]] = dict()
    for i in range(len(tour)):
        if i == 0:
            other_i1 = 1
            other_i2 = -1
        else:
            other_i1 = i+1
            other_i2 = i-1
        
        d[tour[i]] = (tour[other_i1], tour[other_i2])
    
    return d

def random_tour(n: int) -> list[int]:
    randoms = [(random.random(), i) for i in range(n)]
    randoms.sort(key=lambda l: l[0])
    tour = list(map(lambda r: r[1], randoms))
    return tour

def lin_kernighan(n: int, costs: list[list[int]]):
    tour = random_tour(n)
    
    untried_t0 = tour.copy()

    while True:
        i = 0
        t0 = untried_t0.pop()
