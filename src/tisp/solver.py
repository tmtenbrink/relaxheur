import argparse
from tisp.branch_bound import do_branch_and_bound
from tisp.parser import get_inst


def cli():
    default_inst = "tsp/gr48.dat"

    parser = argparse.ArgumentParser(description="TSP solver.")

    nargs = 1 if default_inst is None else "?"

    parser.add_argument(
        "inst_path",
        help="Path to instance.",
        nargs=nargs,
        default=default_inst,
    )

    args = parser.parse_args()

    print("Loading instance...")
    inst, inst_name = get_inst(args.inst_path)

    print(f"Solving inst: {inst_name}")
    print(do_branch_and_bound(inst))


def run():
    cli()


if __name__ == "__main__":
    run()
