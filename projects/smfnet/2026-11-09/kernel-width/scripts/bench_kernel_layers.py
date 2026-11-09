"""Per-layer cost of U-Net convs vs kernel size, at the U-Net levels (8,192
reads), forward + backward (dgrad + wgrad), bf16, channels-last conv2d as
the model runs them; plus the whole conv block (fused per-read norm -> GELU
-> conv) to show the k-independent part. Roofline ideal per conv pass set:
max(memory time at 3.35 TB/s for ~6 tensor passes, compute time at 990
TFLOP/s for 3 x 2 N L C^2 k FLOPs)."""
import sys
import torch
from smf_net.channels_last import conv_channels_last
from smf_net.fused_norm import read_norm_film_gelu
torch.backends.cudnn.benchmark = True
N = 8192
models = {
    "1 bp U-Net": [(1024, 48), (512, 64), (256, 80), (128, 96), (64, 112), (32, 128), (16, 144)],
    "patched U-Net (4 bp)": [(256, 128), (128, 160), (64, 192), (32, 224), (16, 256)],
}
KS = [3, 5, 9, 15, 25]
def bench(fn, n=20):
    for _ in range(3): fn()
    torch.cuda.synchronize(); s, e = torch.cuda.Event(True), torch.cuda.Event(True); s.record()
    for _ in range(n): fn()
    e.record(); torch.cuda.synchronize(); return s.elapsed_time(e) / n
for model, levels in models.items():
    print(f"\n{model}: conv fwd+bwd ms (roofline ideal ms) per kernel size; block = fused norm + conv")
    print(f"{'L':>5s} {'C':>4s} " + " ".join(f"{'k=' + str(k):>15s}" for k in KS) + f" {'block k=5':>10s} {'block k=9':>10s}")
    for L, C in levels:
        x = torch.randn(N, L, C, device="cuda", dtype=torch.bfloat16, requires_grad=True)
        dy = torch.randn_like(x)
        row = []
        for k in KS:
            w = (torch.randn(C, C, k, device="cuda") * 0.02).to(torch.bfloat16).requires_grad_()
            t = bench(lambda: torch.autograd.grad(conv_channels_last(x, w, None), [x, w], dy))
            mem = 6 * x.numel() * 2 / 3.35e12 * 1e3
            flop = 3 * 2 * N * L * C * C * k / 990e12 * 1e3
            row.append(f"{t:6.2f} ({max(mem, flop):5.2f})")
        blocks = []
        g = torch.ones(C, device="cuda", requires_grad=True); b = torch.zeros(C, device="cuda", requires_grad=True)
        for k in (5, 9):
            w = (torch.randn(C, C, k, device="cuda") * 0.02).to(torch.bfloat16).requires_grad_()
            def step():
                h = read_norm_film_gelu(x, g, b, None, 32)
                torch.autograd.grad(conv_channels_last(h, w, None), [x, w, g, b], dy)
            blocks.append(f"{bench(step):10.2f}")
        print(f"{L:5d} {C:4d} " + " ".join(f"{r:>15s}" for r in row) + " " + " ".join(blocks), flush=True)
        del x, dy; torch.cuda.empty_cache()
