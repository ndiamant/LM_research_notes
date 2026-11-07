"""Fused kernel vs compiled PyTorch for per-read norm -> affine -> FiLM -> GELU,
channels-last, forward+backward, at the U-Net's levels (256 regions x 32 reads)."""
import torch
from smf_net.fused_norm import read_norm_film_gelu, reference
def bench(fn, n=20):
    for _ in range(3): fn()
    torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
    s, e = torch.cuda.Event(True), torch.cuda.Event(True); s.record()
    for _ in range(n): fn()
    e.record(); torch.cuda.synchronize(); return s.elapsed_time(e) / n, torch.cuda.max_memory_allocated() / 2**30
B, R = 256, 32; N = B * R
compiled_ref = torch.compile(reference, dynamic=False)
for L, C in [(1024, 48), (512, 64), (256, 80), (128, 96), (64, 112)]:
    x = torch.randn(N, L, C, device="cuda", dtype=torch.bfloat16, requires_grad=True)
    g = torch.ones(C, device="cuda", requires_grad=True); b = torch.zeros(C, device="cuda", requires_grad=True)
    st = (torch.randn(B, L, 2 * C, device="cuda") * 0.1).to(torch.bfloat16).requires_grad_()
    dh = torch.randn_like(x)
    ideal = 5 * x.numel() * 2 / 3.35e12 * 1e3
    res = {}
    for name, fn in [("compiled PyTorch", compiled_ref), ("fused kernel", read_norm_film_gelu)]:
        def step():
            torch.autograd.grad(fn(x, g, b, st, R), [x, g, b, st], dh)
        res[name] = bench(step)
    (tr, mr), (tf, mf) = res["compiled PyTorch"], res["fused kernel"]
    print(f"  L={L:4d} C={C:3d}: compiled {tr:6.2f} ms ({mr:4.1f} GiB peak)  fused {tf:6.2f} ms ({mf:4.1f} GiB)  "
          f"{tr/tf:4.2f}x faster; fused at {ideal/tf:4.0%} of the bandwidth ideal ({ideal:.2f} ms)")
    del x, st, dh; torch.cuda.empty_cache(); torch._dynamo.reset()
