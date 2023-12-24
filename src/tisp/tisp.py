from tisp.branch_bound import do_branch_and_bound
from tisp.parser import get_inst


def solve_from_path(path_name: str):
    inst, inst_name = get_inst(path_name)

    print(f"Solving inst: {inst_name}")
    print(do_branch_and_bound(inst))
