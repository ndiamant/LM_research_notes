"""Consensus motifs for the HMM repo's GM12878 TF panel, for scripts/marginalize.py.

For each PWM in smf_hmm/evaluation/tf_chip/references (the matrices of
configs/gm12878_train.yaml), the argmax base per position, with end positions of
information content < 0.5 bits trimmed; control: of the first 2,000 matched_controls() shuffles (same
assay site counts), the one sharing the shortest substring with the motif in either
orientation, then the most distant (the first shuffle is often a single swap away).
"""
import numpy as np

from itertools import islice

from smf_net.evaluation.marginalization import matched_controls, reverse_complement

REF = "/home/users/diamant/repos/smf_hmm/evaluation/tf_chip/references"
TFS = {"CTCF": "MA0139.2", "NRF1": "MA0506.1", "USF2": "MA0526.3", "CREM": "MA0609.3",
       "RFX5": "MA0510.1", "EBF1": "MA0154.4", "PAX5": "MA0014.2", "SPI1": "MA0080.5", "YY1": "MA0095.2"}
print("name\tsequence\tcontrol_sequence\tjaspar_id\tgroup")
for name, jid in TFS.items():
    rows, on = [], False
    for line in open(f"{REF}/{jid}.meme"):
        if line.startswith("letter-probability"):
            on = True
            continue
        if on:
            vals = line.split()
            if len(vals) != 4:
                if rows:
                    break
                continue
            rows.append([float(v) for v in vals])
    p = np.array(rows)
    ic = 2 + (p * np.log2(np.clip(p, 1e-9, 1))).sum(1)
    keep = np.flatnonzero(ic >= 0.5)
    p = p[keep[0] : keep[-1] + 1]
    seq = "".join("ACGT"[i] for i in p.argmax(1))
    targets = (seq, reverse_complement(seq))

    def shared(c):  # longest substring shared with the motif in either orientation
        return max(k for k in range(len(c) + 1) if k == 0 or any(
            c[i:i + k] in t for i in range(len(c) - k + 1) for t in targets))

    def score(c):
        return (shared(c), -min(sum(a != b for a, b in zip(c, t)) for t in targets))

    control = min(islice(matched_controls(seq), 2000), key=score)
    print(f"{name}\t{seq}\t{control}\t{jid}\tgm12878_hmm")
