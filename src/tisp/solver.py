import argparse
from tisp.branch_bound.branch_bound import do_branch_and_bound
from tisp.lp.model import lp_initialize
from tisp.lp.optimize import lp_optimize_optimal
from tisp.parser import get_inst
from tisp.types import Formulation


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

    # lp = lp_initialize(Formulation.CUTTING_PLANE, inst)

    print(f"Solving inst: {inst_name}")
    print(do_branch_and_bound(inst))


def run():
    cli()


if __name__ == "__main__":
    run()
