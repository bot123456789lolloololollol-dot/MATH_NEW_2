"""Load a champion program from npz and expose a callable `bins(sizes, cap)`."""
import numpy as np
from .gp import pack_gp, program_to_str, ProgSpec


def load_champion(path):
    z = np.load(path)
    ops = z["ops"].astype(np.int32)
    consts = z["consts"].astype(np.float64)
    prog = ProgSpec(ops, consts)

    def bins(sizes, cap):
        s = np.asarray(sizes)
        if s.dtype != np.int32:
            s = s.astype(np.int32)
        sd = np.sort(s)[::-1].copy()
        return int(pack_gp(sd, int(cap), prog.ops, prog.consts)[0])

    bins.prog = prog
    bins.str = program_to_str(prog)
    return bins
