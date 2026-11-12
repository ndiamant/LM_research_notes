"""Whole-molecule m6A rate with read IDs, from the GM12878 fiberseq calls TSVs.

Same counts as SMFNet notes' read_rates.py (accessible '1' / (accessible '1' + protected '2')
over the full aligned span), plus read_id, so H5 rows can be joined to molecule rates.
Usage: molecule_rates_named.py OUT_DIR CHROM...
"""
import subprocess, sys
from multiprocessing import Pool
from pathlib import Path
import pandas as pd

CALLS = Path("/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/calls")
OUT = Path(sys.argv[1])

def one(chrom):
    proc = subprocess.Popen(["/bin/zstd", "-dc", "--long=27", str(CALLS / f"{chrom}.tsv.zst")],
                            stdout=subprocess.PIPE, bufsize=1 << 22)
    rows = []
    for line in proc.stdout:
        if line.startswith(b"#"):
            continue
        f = line.rstrip(b"\n").split(b"\t")
        calls = f[6]
        rows.append((f[3].decode(), int(f[1]), int(f[2]), calls.count(b"1"), calls.count(b"2")))
    if proc.wait() != 0:
        raise RuntimeError(f"zstd failed on {chrom}")
    t = pd.DataFrame(rows, columns=["read_id", "start", "end", "acc", "pro"])
    if t.read_id.duplicated().any():
        raise ValueError(f"duplicate read IDs on {chrom}")
    t.to_csv(OUT / f"{chrom}.tsv.gz", sep="\t", index=False)
    return chrom, len(t)

if __name__ == "__main__":
    with Pool(len(sys.argv) - 2) as pool:
        for chrom, n in pool.imap_unordered(one, sys.argv[2:]):
            print(chrom, n, flush=True)
