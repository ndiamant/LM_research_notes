---
title: Do SMFNet models memorize training regions' DNA? Train vs val loss by mask rate, with own, wrong and mutated DNA
date: 2026-11-08
project: smfnet
agent: Claude Code (claude-opus-5-5)
status: draft
sources:
  - memorization/scripts/memorization_gap.py, run_all.sh (Sherlock job 51525012)
  - memorization/scripts/summary.txt, gap_control.txt
tags: [overfitting, memorization, generalization, alphagenome, mask-rate, evaluation]
---

# Summary

The question: do models with a trained DNA encoder (the original SMFNet,
fine-tuned AlphaGenome, the ATAC model) memorize the DNA of training
regions, while frozen-AlphaGenome models generalize? The training curves
suggested it, because the original SMFNet reaches nearly the lowest train
loss but not the lowest val loss.

Five models were scored identically on 4,000 train and 4,000 val regions
per h5 (active and background), using:
- fixed masks at 11 mask rates;
- the region's own DNA;
- another region's DNA ("rolled");
- the region's DNA with 2% of bases substituted.

**Findings:**

1. **A shared floor.** Every model's val − train gap with the wrong DNA is
   the same, ~0.0035–0.0047 per hidden call at mask rates up to 0.55. It
   doesn't depend on how DNA is encoded, so it is not DNA memorization. It is
   either a difference between the train and val regions in how predictable
   reads are from their own calls, or something read-level that every model
   shares (unresolved).
2. **DNA-attributable memorization.** Above that floor (taken as a fixed
   0.004), the gap with the own DNA is clearly model-dependent:

   | model | excess over the floor, mask 0.05 → 0.85 |
   | --- | --- |
   | original SMFNet (CNN DNA encoder) | +0.005 → +0.007 |
   | ATAC model | +0.004 → +0.005 |
   | fine-tuned AlphaGenome | +0.002 → +0.005 |
   | frozen AlphaGenome | ≤ +0.0013 |
   | frozen AlphaGenome + attention U-Net (Muon+EMA) | ≤ +0.0016 |

   So the models whose DNA encoder is trained memorize their training regions'
   DNA, and freezing AlphaGenome removes most of it.
3. **The hypothesized mask-rate trend holds only weakly.**
   - For the original SMFNet, the gap grows from 0.0084 (mask 0.25) to
     0.0113 (0.85).
   - For fine-tuned AlphaGenome it grows from 0.0057 to 0.0091.
   - For the frozen-AlphaGenome models it does not grow, and it shrinks
     toward full masking (0.0005 for frozen AlphaGenome at mask 1.0).
   - For the original SMFNet, the extra DNA benefit on train grows with the
     mask rate, from 0.005 (0.05) to 0.019 (mask 1.0). The DNA benefit here
     is the loss with another region's DNA minus the loss with the own DNA.
4. **Where the val gains come from.** All five models have nearly the same
   val loss at low mask rates (0.193–0.196 at 0.05). The frozen-AlphaGenome
   models' val advantage is at high mask rates, where DNA matters most:
   0.5180 vs 0.5207 at 0.95. There, the original SMFNet's lower train loss
   didn't carry over.
5. **Full masking against the region mean.** At mask 1.0 the prediction is
   DNA-only. Predicting each call from the mean of the region's other reads
   gives 0.6007 (train) and 0.6002 (val), so the two splits are equally hard
   at the level of marginals.
   - The original SMFNet reaches 0.5940 on train, better than that region
     mean, but 0.6033 on val.
   - Frozen AlphaGenome gives 0.5997 on train and 0.6003 on val.

   Reading: the original model learned the training regions' profiles from
   their own reads (in-sample), and that doesn't transfer.
6. **Brittleness to mutated DNA.** The loss increase from mutating 2% of
   bases, train/val ratio at mask 0.05 → 1.0:

   | model | ratio |
   | --- | --- |
   | original SMFNet | 1.2× → 1.6× |
   | ATAC | 1.1× → 1.6× |
   | fine-tuned AlphaGenome | 1.0× → 1.4× |
   | attention U-Net (Muon+EMA) | 1.0× → 1.3× |
   | frozen AlphaGenome | 1.0× → 1.16× |

   The trained-encoder models are more sensitive to small changes in exactly
   the sequences they trained on.

# Key Points

- **What memorization predicts.** If a model memorizes regions through their
  DNA, it shows:
  - a train − val gap that is larger with the own DNA than with another
    region's DNA;
  - a larger DNA benefit on train than on val;
  - more brittleness to mutations on train;
  - train loss below what the region's other reads can predict, with val
    above it.

  The original SMFNet shows all four; frozen-AlphaGenome models show little
  of any.
- **It isn't the main driver of the remaining val gap at low masking.** At
  low mask rates all models are equal on val. The shared ~0.004 floor is not
  DNA memorization.
- **The ATAC model has the best sampled profile metrics, yet it memorizes
  about as much as the original SMFNet.** Its DNA encoder is the trained CNN.
  Inference: a frozen-AlphaGenome + ATAC model might combine ATAC's gain with
  AlphaGenome's generalization (untested).

# Details

## Method (`memorization_gap.py`)

- **Region subset.** For one h5 and split, a fixed 4,000-region subset
  (`max_regions_per_h5_per_split`, the loader's seeded choice).
- **Same scoring on both splits.** Each region is set up as validation is
  (`_setup_split(split, train=False)`): 32 reads chosen by region index, the
  region's middle (no jitter), no reverse complement. Masks are a uniform
  draw per call seeded by batch index, hidden where draw < mask rate; at
  rate 1.0 all valid calls are hidden.
- **Measure.** Cross-entropy per hidden call, pooled over CpG, GpC and GCG
  sites and over active plus background regions.
- **DNA conditions:**
  - own;
  - rolled: the DNA (and ATAC) of the next region in the batch;
  - mutated: 2% of the known bases changed to a random other base, and ATAC
    kept for the ATAC model.
- **Region-mean reference.** For each call, the mean of the same site over
  the region's other 99 of 100 reads (leave-one-out, Jeffreys prior).
  - It is not a strict bound: the per-site estimate is noisy, so a model
    that smooths across neighbouring sites can beat it.
  - It is the same reference on both splits.
- **Models:**
  - original SMFNet: `SMFNet_full_run_2026-08-11` epoch 377;
  - frozen and fully fine-tuned AlphaGenome (`ONT_mESC_alphagenome_frozen`,
    `_ftall_run2`): epoch 377;
  - attention U-Net Muon+EMA: epoch 377;
  - ATAC model (`ONT_mESC_ATAC_TKO_pooled`, nd_atac code): epoch 377.
- **Caveats.**
  - Every train region was trained on with all its reads
    (`train_val_read_fraction=1`), so train scoring is in-sample for reads
    too.
  - Sums are stored pooled, not per region, so there are no standard errors.
    The pattern is consistent across 11 mask rates and 5 models.

## Train − val gap, own DNA / rolled DNA (val − train, nats per hidden call)

| mask | original SMFNet | frozen AG | fine-tuned AG | attn U-Net Muon+EMA | ATAC |
| --- | --- | --- | --- | --- | --- |
| 0.05 | +0.0091 / +0.0043 | +0.0053 / +0.0046 | +0.0064 / +0.0047 | +0.0049 / +0.0047 | +0.0084 / +0.0046 |
| 0.25 | +0.0084 / +0.0035 | +0.0044 / +0.0035 | +0.0057 / +0.0035 | +0.0042 / +0.0035 | +0.0078 / +0.0038 |
| 0.45 | +0.0092 / +0.0040 | +0.0049 / +0.0042 | +0.0067 / +0.0041 | +0.0049 / +0.0042 | +0.0086 / +0.0046 |
| 0.65 | +0.0101 / +0.0025 | +0.0050 / +0.0030 | +0.0078 / +0.0026 | +0.0051 / +0.0030 | +0.0091 / +0.0033 |
| 0.85 | +0.0113 / −0.0016 | +0.0045 / −0.0005 | +0.0091 / −0.0015 | +0.0050 / −0.0004 | +0.0093 / −0.0000 |
| 1.00 | +0.0093 / −0.0101 | +0.0005 / −0.0074 | +0.0065 / −0.0095 | +0.0016 / −0.0076 | +0.0058 / −0.0072 |

**Why the rolled-DNA gap turns negative at high mask rates** (inferred): with
little read context, a model confident about the wrong region's profile is
penalized more on train regions, whose DNA it has memorized. So the
rolled-DNA control is clean only at low to mid mask rates.

The full tables, with train and val losses, DNA benefit and mutation cost per
mask rate, are in [summary.txt](memorization/scripts/summary.txt).

## Commands

```bash
sbatch $SCRATCH/smfnet_memorization/run_all.sbatch   # btrippe, 1 GPU, 6 h limit; ran 1 h 38 min
module load devel pixi/0.53.0; export PIXI_CACHE_DIR=$SCRATCH/pixi-cache
cd $SCRATCH/smfnet_memorization
PIXI="pixi run --frozen --manifest-path $HOME/repos/SMFNet_public/pyproject.toml"
$PIXI python summarize.py; $PIXI python gap_control.py
```

- **Job history.** The first submission (51524753) failed in 12 s: bash
  4.2 with `set -u` rejected an empty array expansion. It was fixed as
  `${ck[@]+"${ck[@]}"}` and resubmitted as 51525012.
- **Code.** Non-ATAC models were run with
  `PYTHONPATH=~/repos/SMFNet_patch/src`, the ATAC model with
  `$SCRATCH/src/SMFNet_atac/src`.
- **Results:** `$SCRATCH/smfnet_memorization/results/<model>_<ac|bg>_<train|val>.json`.

# Related Notes

- [Frozen AlphaGenome encoder](../2026-10-29/frozen-alphagenome-encoder.md):
  the original, frozen, fine-tuned and ATAC models compared here, and their
  val metrics.
- [GCH overfitting, regularization, EMA, and sampling-length experiments](../2026-10-07/regularization-ema-sampling-experiments.md):
  earlier overfitting work on another dataset.
- [Patched U-Net](../2026-11-07/patched-unet.md): the attention U-Net
  Muon+EMA model scored here, and the patched run still training.

# Open Questions

- **What is the shared ~0.004 floor?** A train/val region difference, or
  read-level memorization that every model does (train reads were all seen
  in training)? A model trained with `train_val_read_fraction < 1`, holding
  out reads within training regions, would separate read from region
  memorization.
- **Can a trained DNA encoder be kept from memorizing?** For example, DNA
  mutation or dropout augmentation, or freezing.
- **Would frozen AlphaGenome + ATAC keep ATAC's sampled-profile gains
  without the memorization?**
- **The patched Muon+EMA run** should be added to this comparison when it
  finishes.

# Sources

- [memorization_gap.py](memorization/scripts/memorization_gap.py),
  [run_all.sh](memorization/scripts/run_all.sh),
  [summarize.py](memorization/scripts/summarize.py),
  [gap_control.py](memorization/scripts/gap_control.py); outputs
  [summary.txt](memorization/scripts/summary.txt) and
  [gap_control.txt](memorization/scripts/gap_control.txt).
