"""Per-molecule m6A rate (accessible / (accessible + protected)) from fiberseq calls TSVs."""
import subprocess, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np

CALLS = Path("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/calls")
OUT = Path(sys.argv[1])

def one(path):
    proc = subprocess.Popen(["zstd", "-dc", "--long=27", str(path)], stdout=subprocess.PIPE, bufsize=1 << 22)
    acc, pro, span = [], [], []
    for line in proc.stdout:
        if line.startswith(b"#"):
            continue
        calls = line.rstrip(b"\n").rsplit(b"\t", 1)[1]
        acc.append(calls.count(b"1")); pro.append(calls.count(b"2")); span.append(len(calls))
    if proc.wait() != 0:
        raise RuntimeError(f"zstd failed on {path}")
    np.savez(OUT / f"{path.name.split('.')[0]}.npz", acc=np.array(acc, np.int32),
             pro=np.array(pro, np.int32), span=np.array(span, np.int32))
    return path.name, len(acc)

if __name__ == "__main__":
    files = sorted(CALLS.glob("*.tsv.zst"), key=lambda p: -p.stat().st_size)
    if len(sys.argv) > 2:
        files = [CALLS / f"{c}.tsv.zst" for c in sys.argv[2:]]
    with Pool(8) as pool:
        for name, n in pool.imap_unordered(one, files):
            print(name, n, flush=True)
