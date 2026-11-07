"""Attention U-Net (frozen AG encoder, per-read norm, 2 attention layers +
cross-attention): unfused vs unet_fused_norm. Equivalence on a small batch
(same weights), then compiled training-step time and memory at full batch."""
import sys, time
from types import SimpleNamespace
import torch
from smf_net.cnn import CNNConfig
from smf_net.model import SMFNet, SMFNetConfig

torch.manual_seed(0)
torch.backends.cudnn.benchmark = "--cudnn-benchmark" in sys.argv
R, L = 32, 2048
base = dict(dna_encoder_type="alphagenome", smf_body_type="unet", unet_norm="read",
            unet_attention_layers=2, unet_cross_attention=True)
def build(fused, cl=False):
    torch.manual_seed(0)
    m = SMFNet(SMFNetConfig(cnn_config=CNNConfig(**base, unet_fused_norm=fused, unet_channels_last=cl, output_conv_matmul=cl))).cuda()
    return m
ref, fused = build(False), build(True, cl=True)
fused.load_state_dict(ref.state_dict())
# non-trivial FiLM and affine parameters, so they're exercised
with torch.no_grad():
    for m in (ref, fused):
        torch.manual_seed(1)
        for name, p in m.named_parameters():
            if "cond_proj.1.weight" in name or name.endswith("norm.weight"):
                p.add_(torch.randn_like(p) * 0.05)
def data(b):
    g = torch.Generator().manual_seed(2)
    dna = torch.nn.functional.one_hot(torch.randint(4, (b, L), generator=g), 4).transpose(1, 2).float().cuda()
    smf = torch.nn.functional.one_hot(torch.randint(0, 3, (b, R, L // 2), generator=g), 3).permute(0, 1, 3, 2).float().cuda()
    return smf, dna
smf, dna = data(8)
outs, grads = [], []
for m in (ref, fused):
    m.zero_grad()
    with torch.autocast("cuda", dtype=torch.bfloat16):
        out = m.model(smf, dna).float()
    out.square().mean().backward()
    outs.append(out.detach())
    grads.append({n: p.grad.detach().float().clone() for n, p in m.named_parameters() if p.grad is not None})
rel = lambda a, b: ((a - b).norm() / b.norm().clamp_min(1e-12)).item()
print(f"output rel diff {rel(outs[1], outs[0]):.2e}")
worst = sorted(((rel(grads[1][n], g), n) for n, g in grads[0].items()), reverse=True)[:3]
print("largest gradient rel diffs:", [(f"{e:.1e}", n) for e, n in worst])
print("median gradient rel diff:", f"{torch.tensor([rel(grads[1][n], g) for n, g in grads[0].items()]).median():.1e}")

smf, dna = data(256)
for name, m in (("unfused", ref), ("fused + channels-last + matmul output conv", fused)):
    m._trainer = SimpleNamespace(max_epochs=400, estimated_stepping_batches=182000)
    opt = m.configure_optimizers()["optimizer"]
    m.model.compile()
    torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
    def step():
        with torch.autocast("cuda", dtype=torch.bfloat16):
            loss = m.model(smf, dna).float().square().mean()
        loss.backward(); opt.step(); opt.zero_grad(set_to_none=True)
    for _ in range(5): step()
    torch.cuda.synchronize(); t = time.time()
    for _ in range(15): step()
    torch.cuda.synchronize()
    print(f"{name}: {(time.time() - t) / 15 * 1000:.0f} ms per step, peak {torch.cuda.max_memory_allocated() / 2**30:.1f} GiB (cudnn.benchmark={torch.backends.cudnn.benchmark})")
