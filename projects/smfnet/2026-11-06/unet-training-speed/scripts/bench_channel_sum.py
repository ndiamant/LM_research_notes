"""Column sum of (rows, C) bf16 -> fp32 (the channels-last bias gradient):
native torch sum vs a Triton two-stage kernel, at the U-Net level sizes."""
import torch, triton, triton.language as tl

@triton.jit
def _partial(X, P, rows, C, ROWS_PER: tl.constexpr, BLOCK_R: tl.constexpr, BLOCK_C: tl.constexpr):
    pid = tl.program_id(0)
    cols = tl.arange(0, BLOCK_C)
    cmask = cols < C
    acc = tl.zeros([BLOCK_R, BLOCK_C], dtype=tl.float32)
    for start in range(pid * ROWS_PER, (pid + 1) * ROWS_PER, BLOCK_R):
        r = start + tl.arange(0, BLOCK_R)
        m = (r[:, None] < rows) & cmask[None, :]
        acc += tl.load(X + r[:, None] * C + cols[None, :], mask=m, other=0.0).to(tl.float32)
    tl.store(P + pid * C + cols, tl.sum(acc, axis=0), mask=cmask)

def triton_sum(x, rows_per=4096, block_r=64, warps=4):
    x2 = x.reshape(-1, x.shape[-1]); rows, C = x2.shape
    progs = triton.cdiv(rows, rows_per)
    p = torch.empty(progs, C, device=x.device, dtype=torch.float32)
    _partial[(progs,)](x2, p, rows, C, ROWS_PER=rows_per, BLOCK_R=block_r, BLOCK_C=triton.next_power_of_2(C), num_warps=warps)
    return p.sum(0)

def bench(f, n=30):
    for _ in range(3): f()
    torch.cuda.synchronize(); s, e = torch.cuda.Event(True), torch.cuda.Event(True); s.record()
    for _ in range(n): f()
    e.record(); torch.cuda.synchronize(); return s.elapsed_time(e) / n

N = 8192
for L, C in [(1024, 48), (512, 64), (256, 80), (128, 96), (64, 112), (32, 128), (16, 144), (1024, 1)]:
    x = torch.randn(N, L, C, device="cuda", dtype=torch.bfloat16)
    ideal = x.numel() * 2 / 3.35e12 * 1e3
    ref = x.reshape(-1, C).sum(0, dtype=torch.float32)
    out = [f"native {bench(lambda: x.reshape(-1, C).sum(0, dtype=torch.float32)):.3f}"]
    for rp, br, w in [(2048, 64, 4), (4096, 64, 4), (8192, 128, 8), (16384, 128, 8)]:
        got = triton_sum(x, rp, br, w)
        assert torch.allclose(got, ref, rtol=1e-3, atol=1e-1), (got - ref).abs().max()
        out.append(f"{rp}/{br}/{w} {bench(lambda: triton_sum(x, rp, br, w)):.3f}")
    print(f"L={L:4d} C={C:3d} ideal {ideal:.3f} | " + "  ".join(out) + " ms", flush=True)
