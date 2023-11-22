import argparse
from pathlib import Path

def parse_line(ln: str):
    return list(map(lambda i: int(i), ln.rstrip().split(" ")))


def get_inst_path():
    parser = argparse.ArgumentParser(description="Process some integers.")

    parser.add_argument("inst_name", help="Filename of instance")

    args = parser.parse_args()

    return Path(args.inst_name)


def parse_as_adj_matrix(inst_path: Path):
    with open(inst_path, "r") as f:
        lines = f.readlines()

    line_0 = parse_line(lines[0])
    n = line_0[0]
    adj_matrix = list(map(lambda ln: parse_line(ln), lines[1:]))

    return n, adj_matrix


def run():
    # inst_path = get_inst_path()
    inst_path = Path('tsp/gr48.dat')

    n, graph_l = parse_as_adj_matrix(inst_path)

    print(graph_l)
    print(n)


if __name__ == "__main__":
    run()
