"""One-output-channel conv (C -> 1, k=75) as tensor-core matmuls: for each block
of BL output positions, y[n, l0:l0+BL] = h_window[n] (flattened (BL+k-1) x C)
@ T, with T[(r, c), i] = w[c, r - i] for 0 <= r - i < k (a Toeplitz matrix).
Channels-last h (N, L, C)."""
import torch, torch.nn.functional as F
from smf_net.cnn import conv_to_one_channel
torch.backends.cudnn.benchmark = True
def toeplitz(w, block):
    # T[(r, c), i] = w[c, r - i]: (block + k - 1, C, block) -> ((block + k - 1) * C, block)
    _, C, k = w.shape
    rows = block + k - 1
    r = torch.arange(rows, device=w.device)[:, None]; i = torch.arange(block, device=w.device)[None, :]
    j = r - i                                               # tap index
    valid = (j >= 0) & (j < k)
    T = w[0][:, j.clamp(0, k - 1)] * valid                  # (C, rows, block)
    return T.permute(1, 0, 2).reshape(rows * C, block)
def conv_toeplitz(h, w, b, block):
    N, L, C = h.shape; k = w.shape[-1]; pad = k // 2
    hp = F.pad(h, (0, 0, pad, pad))                          # (N, L + k - 1, C)
    T = toeplitz(w, block).to(h.dtype)
    rows = block + k - 1
    outs = [hp[:, l0:l0 + rows].reshape(N, rows * C) @ T for l0 in range(0, L, block)]
    return torch.cat(outs, dim=1) + b.to(h.dtype)
def bench(fn, n=10):
    for _ in range(3): fn()
    torch.cuda.synchronize(); s, e = torch.cuda.Event(True), torch.cuda.Event(True); s.record()
    for _ in range(n): fn()
    e.record(); torch.cuda.synchronize(); return s.elapsed_time(e) / n
torch.backends.cuda.matmul.allow_tf32 = False; torch.backends.cudnn.allow_tf32 = False
h = torch.randn(4, 1024, 48, device="cuda"); w = torch.randn(1, 48, 75, device="cuda"); b = torch.randn(1, device="cuda")
ref = F.conv1d(h.transpose(1, 2), w, b, padding=37)[:, 0]
for block in (64, 128):
    print(f"  block {block}: rel diff vs conv1d (fp32) {((conv_toeplitz(h, w, b, block) - ref).norm() / ref.norm()).item():.1e}")
for N, L, C in [(8192, 1024, 48), (8192, 1024, 128)]:
    h = torch.randn(N, L, C, device="cuda", dtype=torch.bfloat16, requires_grad=True)
    w = (torch.randn(1, C, 75, device="cuda") * 0.02).requires_grad_(); b = torch.zeros(1, device="cuda", requires_grad=True)
    res = {}
    for name, fn in [("matmul + shift", lambda h, w, b: conv_to_one_channel(h, w, b, True)[:, 0])] + \
                    [(f"Toeplitz matmul, block {bl}", (lambda bl: lambda h, w, b: conv_toeplitz(h, w, b, bl))(bl)) for bl in (128, 256, 512)]:
        cf = torch.compile(fn)
        def step():
            with torch.autocast("cuda", dtype=torch.bfloat16): y = cf(h, w, b)
            y.float().sum().backward()
        res[name] = bench(step); torch._dynamo.reset()
    print(f"  N={N} L={L} C={C}: " + "   ".join(f"{k} {v:.2f} ms" for k, v in res.items()))
