"""Per-read and per-read-per-site tables for the read m6A-rate evaluation, one val chromosome.

Usage: build_read_tables.py CHROM OUT_DIR

Reads: every read of every benchmark window, with its whole-molecule m6A rate (calls TSV,
joined through active.h5 read_id), its rate over the window's informative A/T calls (±3 bp
context efficiency >= 0.25 in the middle openness bin of the bias table), and each model's
mean posterior occupancy per state group over the 512 bp window.

Read-sites: every read at every IDR-labelled unique motif site (analysis_idr), with the read
class at the motif (footprint / open / closed / other, as in read_classes.py), whether it has
any protected call in the motif (the analysis's protected-read definition), and each model's
mean P(TF) over the motif.

Also accumulated per molecule-rate bin (BINS): the m6A autocorrelation along reads at lags
0..MAX_LAG, and m6A profiles (strand-oriented, +-PROFILE bp) around IDR-bound and unbound
CTCF and NRF1 motifs.
"""

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, "/scratch/users/diamant/fiber_eval_work/tf_screen")
from genome_pilot import kmer_codes  # noqa: E402

CHROM, OUT = sys.argv[1], Path(sys.argv[2])
ROOT = Path(f"/scratch/users/diamant/gm12878_validation/covariation_windows_v2/{CHROM}-512-valid90-reads10-flank128")
H5 = "/scratch/groups/btrippe/ndiamant/fiberseq/GM12878/h5_ont_polarity/active.h5"
RATES = f"/scratch/users/diamant/fiber_read_m6a/named/{CHROM}.tsv.gz"
MODELS = ["phase_type", "ag_kl10", "ag_ol001", "ag_ol001_s05", "ag_kl10_s0"]
FLANK = (15, 40)
BINS = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 1.01])
MAX_LAG = 300
PROFILE = 125
PROFILE_TFS = ("ctcf", "nrf1")
table = np.load("/scratch/users/diamant/fiber_eval_work/tf_screen/bias_table_chr8_chr20.npz")
EFF = (1 - table["exp_protection"][5]) / table["bin_rate"][5]


def read_table(group):
    return pd.DataFrame({c: (group[c].asstr()[:] if h5py.check_string_dtype(group[c].dtype) else group[c][:])
                         for c in json.loads(group.attrs["columns"])})


def tf_index(names):
    return [i for i, n in enumerate(names) if n == "TF" or n.startswith("TF_")]


rates = pd.read_csv(RATES, sep="\t", usecols=["read_id", "acc", "pro"]).set_index("read_id")
sites = pd.read_csv(ROOT / "analysis_idr/unique_motif_summaries.tsv.gz", sep="\t",
                    usecols=["site_key", "tf", "window_id", "start", "end", "center", "strand", "idr_label",
                             "chip_pooled_sum", "accessibility", "window_start"])
sites_by_window = {w: g for w, g in sites.groupby("window_id")}

source = h5py.File(ROOT / "dataset.h5")
meta = read_table(source["window_metadata"]).set_index("window_id")
posteriors = {m: h5py.File(ROOT / f"posteriors/{m}.h5") for m in MODELS}
names = {m: list(f["group_names"].asstr()[:]) for m, f in posteriors.items()}
window_ids = list(posteriors[MODELS[0]]["window_ids"].asstr()[:])
offsets = posteriors[MODELS[0]]["window_offsets"][:]
for m, f in posteriors.items():
    assert list(f["window_ids"].asstr()[:]) == window_ids and np.array_equal(f["window_offsets"][:], offsets), m
h5 = h5py.File(H5)
region_ids = {}

n_bins = len(BINS) - 1
acf_num = np.zeros((n_bins, MAX_LAG + 1))
acf_den = np.zeros((n_bins, MAX_LAG + 1))
prof_m6a = {}
prof_called = {}
read_rows, site_rows = [], []
missing = 0
LIMIT = int(sys.argv[3]) if len(sys.argv) > 3 else None
for wi, wid in enumerate(window_ids[:LIMIT]):
    w = source["windows"][wid]
    dna, pos, obs = w["dna"][:].astype(np.int8), w["positions"][:], w["observations"][:]
    mol = w["molecule_id"].asstr()[:]
    n, L = len(obs), len(dna)
    calls = np.full((n, L), -1, np.int8)
    calls[:, pos] = obs
    codes = kmer_codes(dna)
    informative = (codes >= 0) & (EFF[np.maximum(codes, 0)] >= 0.25)
    inf_calls = np.where(informative[None, :], calls, -1)
    # Whole-molecule rate through active.h5 region:row -> read_id.
    read_ids = []
    for m_id in mol:
        region, row = m_id.rsplit(":", 1)
        if region not in region_ids:
            region_ids[region] = h5["read_id"][region][:].astype(str)
        read_ids.append(region_ids[region][int(row)])
    joined = rates.reindex(read_ids)
    missing += int(joined.acc.isna().sum())
    mol_rate = (joined.acc / (joined.acc + joined.pro)).to_numpy()
    called = inf_calls >= 0
    m6a = inf_calls == 0
    win_called, win_m6a = called.sum(1), m6a.sum(1)
    a, b = offsets[wi], offsets[wi + 1]
    post = {m: posteriors[m]["posterior"][a:b] for m in MODELS}
    row = dict(window_id=wid, molecule_id=mol, read_id=read_ids, mol_rate=mol_rate,
               win_called=win_called, win_m6a=win_m6a)
    for m in MODELS:
        occ = post[m].mean(1)
        for gi, g in enumerate(names[m]):
            row[f"{m}_{g}"] = occ[:, gi]
    read_rows.append(pd.DataFrame(row))

    # Autocorrelation of m6A along each read, standardized by the read's window rate.
    bin_idx = np.digitize(np.nan_to_num(mol_rate, nan=-1), BINS) - 1
    r = np.where(win_called > 0, win_m6a / np.maximum(win_called, 1), np.nan)
    ok = (bin_idx >= 0) & (bin_idx < n_bins) & (win_called >= 20) & (r > 0) & (r < 1)
    if ok.any():
        x = np.where(called[ok], m6a[ok] - r[ok, None], 0.0)
        c = called[ok].astype(float)
        size = 1 << int(np.ceil(np.log2(2 * L)))
        fx, fc = np.fft.rfft(x, size), np.fft.rfft(c, size)
        num = np.fft.irfft(fx * np.conj(fx), size)[:, : MAX_LAG + 1]
        den = np.fft.irfft(fc * np.conj(fc), size)[:, : MAX_LAG + 1] * (r[ok] * (1 - r[ok]))[:, None]
        np.add.at(acf_num, bin_idx[ok], num)
        np.add.at(acf_den, bin_idx[ok], den)

    if wid not in sites_by_window:
        continue
    ws = int(meta.loc[wid, "window_start"])
    for s in sites_by_window[wid].itertuples():
        m0, m1 = s.start - ws, s.end - ws
        Lf = inf_calls[:, max(m0 - FLANK[1], 0):m0 - FLANK[0]]
        Rf = inf_calls[:, m1 + FLANK[0]:m1 + FLANK[1]]
        M = inf_calls[:, m0:m1]
        m_called, m_m6a = (M >= 0).sum(1), (M == 0).sum(1)
        l_m6a, r_m6a = (Lf == 0).sum(1), (Rf == 0).sum(1)
        flanks_open = (l_m6a >= 2) & (r_m6a >= 2)
        cls = np.full(n, "other", object)
        cls[(m_called >= 2) & (m_m6a == 0) & flanks_open] = "footprint"
        cls[(m_m6a >= 1) & flanks_open] = "open"
        cls[(m_called >= 2) & (m_m6a == 0) & (l_m6a == 0) & (r_m6a == 0)
            & ((Lf >= 0).sum(1) >= 2) & ((Rf >= 0).sum(1) >= 2)] = "closed"
        srow = dict(site_key=s.site_key, tf=s.tf, idr_label=s.idr_label, window_id=wid, read_id=read_ids,
                    mol_rate=mol_rate, cls=cls, protected_in_motif=(calls[:, m0:m1] == 1).any(1),
                    motif_called=m_called)
        for m in MODELS:
            srow[m] = post[m][:, m0:m1][..., tf_index(names[m])].sum(-1).mean(1)
        site_rows.append(pd.DataFrame(srow))
        if s.tf in PROFILE_TFS and s.idr_label in ("bound", "unbound"):
            c0 = s.center - ws
            lo, hi = c0 - PROFILE, c0 + PROFILE + 1
            if lo < 0 or hi > L:
                continue
            seg = inf_calls[:, lo:hi]
            if s.strand == "-":
                seg = seg[:, ::-1]
            key = (s.tf, s.idr_label)
            if key not in prof_m6a:
                prof_m6a[key] = np.zeros((n_bins, hi - lo))
                prof_called[key] = np.zeros((n_bins, hi - lo))
            valid = (bin_idx >= 0) & (bin_idx < n_bins)
            np.add.at(prof_m6a[key], bin_idx[valid], (seg[valid] == 0).astype(float))
            np.add.at(prof_called[key], bin_idx[valid], (seg[valid] >= 0).astype(float))
    if wi % 500 == 0:
        print(CHROM, wi, len(window_ids), flush=True)

OUT.mkdir(parents=True, exist_ok=True)
pd.concat(read_rows, ignore_index=True).to_csv(OUT / f"reads_{CHROM}.tsv.gz", sep="\t", index=False)
pd.concat(site_rows, ignore_index=True).to_csv(OUT / f"read_sites_{CHROM}.tsv.gz", sep="\t", index=False)
np.savez(OUT / f"accumulators_{CHROM}.npz", bins=BINS, acf_num=acf_num, acf_den=acf_den,
         **{f"prof_m6a_{t}_{l}": v for (t, l), v in prof_m6a.items()},
         **{f"prof_called_{t}_{l}": v for (t, l), v in prof_called.items()})
print(f"{CHROM}: {sum(len(r) for r in read_rows):,} reads, {missing:,} without a molecule rate")
