"""Forward-only: two-kernel (stats, then apply) vs single-pass fused-norm forward,
channels-last, at the U-Net levels (256 regions x 32 reads); also checks they agree."""
import torch
import smf_net.fused_norm as fn
def bench(f, n=30):
    for _ in range(3): f()
    torch.cuda.synchronize(); s, e = torch.cuda.Event(True), torch.cuda.Event(True); s.record()
    for _ in range(n): f()
    e.record(); torch.cuda.synchronize(); return s.elapsed_time(e) / n
B, R = 256, 32; N = B * R
for L, C in [(1024, 48), (512, 64), (256, 80), (128, 96), (64, 112), (32, 128), (16, 144)]:
    x = torch.randn(N, L, C, device="cuda", dtype=torch.bfloat16)
    g = torch.randn(C, device="cuda"); b = torch.randn(C, device="cuda")
    st = (torch.randn(B, L, 2 * C, device="cuda") * 0.1).to(torch.bfloat16)
    ideal = 2 * x.numel() * 2 / 3.35e12 * 1e3
    row = []
    outs = []
    for single, block, warps in [(False, 0, 0), (True, 1024, 4), (True, 2048, 8), (True, 4096, 8), (True, 4096, 16)]:
        fn.SINGLE_PASS_FORWARD, fn.SINGLE_PASS_BLOCK, fn.SINGLE_PASS_WARPS = single, block or 2048, warps or 8
        with torch.no_grad():
            outs.append(fn.read_norm_film_gelu(x, g, b, st, R))
            row.append(f"{bench(lambda: fn.read_norm_film_gelu(x, g, b, st, R)):.3f}")
    err = max((o.float() - outs[0].float()).abs().max().item() for o in outs[1:])
    print(f"L={L:4d} C={C:3d} ideal {ideal:.3f} | two-kernel {row[0]} | single 1024/4 {row[1]}  2048/8 {row[2]}  4096/8 {row[3]}  4096/16 {row[4]} ms | max diff {err:.3g}", flush=True)
