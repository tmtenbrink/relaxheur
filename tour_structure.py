from random import random
from typing import Iterable


def normalize_tour_repr(tour: list[int]):
    max_size = len(tour)

    bytes_per_n = (max_size.bit_length() + 7) // 8

    i_zero = tour.index(0)
    prev_part = tour[:i_zero]
    normalized_tour = tour[i_zero:] + prev_part
    byte_seq = [i.to_bytes(bytes_per_n, byteorder='big') for i in normalized_tour]
    tour_bytes = b''.join(byte_seq)
    tour_bytes_rev = byte_seq[0] + b''.join(reversed(byte_seq[1:]))

    return min(tour_bytes, tour_bytes_rev)


# def random_tour(n: int) -> list[int]:
#     return shuffle_iter(range(n))


# def shuffle_iter(itrble: Iterable[int]) -> list[int]:
#     """Returns a shuffled copy of the iterable as a list."""
#     randoms = [(random(), i) for i in itrble]
#     randoms.sort(key=lambda l: l[0])
#     shuffled = list(map(lambda r: r[1], randoms))
#     return shuffled

# n = 8

# tours = [shuffle_iter(tour_n) for tour_n in [list(range(n)) for _ in range(20)]]

# for t in tours:
#     rev = normalize_tour_repr(list(reversed(t + [3])))
#     not_rev = normalize_tour_repr(t)
#     print(rev == not_rev)

# print()