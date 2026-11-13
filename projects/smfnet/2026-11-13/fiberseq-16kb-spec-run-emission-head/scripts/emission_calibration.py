"""Is the spec model's fixed epsilon floor miscalibrated? On val items (no molecule filter),
fully masked reads: per b_m quintile and network accessibility q bin, the mean predicted
P(m6A) against the observed m6A fraction."""
import os, sys
import numpy as np, torch
from smf_net.diffusion_schedule import MASK_TOKEN
from smf_net.figwig_data import FigwigDataModuleConfig, FigwigSMFDataModule
from smf_net.model import load_model

ckpt = sys.argv[1]
F_DIR = os.environ["SCRATCH"] + "/fiberseq/GM12878"
cfg = FigwigDataModuleConfig(
    fasta_path="/scratch/groups/btrippe/ndiamant/reference/hg38.fa",
    figwig_dir=F_DIR + "/figwig", regions_bed=F_DIR + "/kept_windows.bed",
    assay_type="fiberseq_m6a", region_length=18498, model_input_length=18432, dna_flank=1024,
    min_read_fraction=0.75, val_chrm=["chr8", "chr20"], batch_size=8, num_workers=4,
    max_regions_per_split=160, efficiency_dir=F_DIR + "/efficiency")
data = FigwigSMFDataModule(cfg); data.setup("validate")
model = load_model(ckpt).cuda().eval()
em = model.emission
print(f"learned bias {em.bias.item():+.3f} (e at median b, mean context {torch.sigmoid(em.bias).item():.3f}), "
      f"efficiency scale {em.efficiency_scale.item():.3f}, epsilon {em.epsilon}")
rows = []
with torch.no_grad():
    for batch in data.val_dataloader():
        dna = batch.dna.cuda().float(); eff = batch.efficiency.cuda()
        x0 = batch.smf.cuda().long()[..., 1024:-1024]
        x_t = torch.full_like(x0, MASK_TOKEN)
        onehot = model._encode_smf_tokens(x_t)
        # the network's accessibility logits, as _predict computes them
        from smf_net.model import one_hot_to_tokens
        tokens = one_hot_to_tokens(dna)
        ctx = model._center_crop_smf(em.site_context(tokens), tokens.shape[-1])
        lab = em.label_logit(ctx, eff)
        sites = model._center_crop_smf(model.cfg and __import__("smf_net.data", fromlist=["x"]).get_valid_mask(tokens, "fiberseq_m6a"), tokens.shape[-1])
        x = torch.cat([onehot, (torch.sigmoid(lab) * sites[:, None])[:, :, None].to(onehot.dtype)], 2)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            z = model.model(x, dna=dna, read_cond=eff[..., None]).float().squeeze(2)
        q = torch.sigmoid(z); e = torch.sigmoid(lab)
        p = q * e + (1 - q) * em.epsilon
        obs = x0 != -1
        b = eff[..., None].expand_as(x0)
        rows.append(torch.stack([b[obs], q[obs], p[obs], (x0[obs] == 0).float()], 1).cpu())
r = torch.cat(rows).numpy()
bq = np.quantile(r[:, 0], [0, .2, .4, .6, .8, 1])
print(f"{len(r):,} sites")
print("b quintile          | q<0.05: share, pred P(m6A), obs  | all sites: pred, obs")
for k in range(5):
    s = (r[:, 0] >= bq[k]) & (r[:, 0] <= bq[k + 1]); lo = s & (r[:, 1] < 0.05)
    print(f"{bq[k]:+.2f} to {bq[k+1]:+.2f}      | {lo.sum()/s.sum():.2f}  {r[lo,2].mean():.4f}  {r[lo,3].mean():.4f} | {r[s,2].mean():.4f} {r[s,3].mean():.4f}")
