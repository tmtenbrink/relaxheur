from pathlib import Path

from tisp.types import Costs
from tisp.resources import res_path


def parse_line(ln: str) -> list[float]:
    return list(map(lambda i: float(i), ln.rstrip().split(" ")))


def get_path(inst_path: str):
    base_path = Path(inst_path)

    if base_path.exists():
        return base_path

    from_res_path = res_path.joinpath(inst_path)

    if from_res_path.exists():
        return from_res_path

    from_res_path_tsp = res_path.joinpath("tsp").joinpath(inst_path)

    if from_res_path_tsp.exists():
        return from_res_path_tsp

    raise FileNotFoundError(
        f"Could not find instance path. Looked at path {inst_path}, {from_res_path} and {from_res_path_tsp}"
    )


def get_inst(inst_path: str) -> tuple[Costs, str]:
    path = get_path(inst_path)
    path_name = path.stem

    lines = read_inst_to_lines(path)

    return parse_as_adj_matrix(lines), path_name


def read_inst_to_lines(inst_path: Path) -> list[str]:
    with open(inst_path, "r") as f:
        lines = f.readlines()

    return lines


def parse_as_adj_matrix(lines: list[str]) -> Costs:
    cost_adj_matrix = list(map(lambda ln: parse_line(ln), lines[1:]))

    return cost_adj_matrix
