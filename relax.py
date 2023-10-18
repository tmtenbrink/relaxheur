import argparse
from pathlib import Path
import numpy as np

def parse_line(ln: str):
    return list(map(lambda i: int(i), ln.rstrip().split(' ')))


def run():
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('inst_name', help='Filename of instance')

    args = parser.parse_args()

    inst_path = Path(args.inst_name)

    with open(inst_path, "r") as f:
        lines = f.readlines()

    line_0 = parse_line(lines[0])
    n = line_0[0]
    all_lines = list(map(lambda ln: parse_line(ln), lines[1:]))
    array = np.array(all_lines)

    print(array)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
