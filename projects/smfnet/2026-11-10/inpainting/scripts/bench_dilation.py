"""Cost of dilated vs undilated convs at the attention U-Net's finest levels
(8,192 reads), forward + backward, bf16, as conv2d on channels-last (N, C, 1,
L) views like the model's convs; plus the fused norm -> conv block for scale."""
import torch
import torch.nn.functional as F
from smf_net.fused_norm import read_norm_film_gelu
torch.backends.cudnn.benchmark = True
N = 8192
def conv(h, w, d):
    k = w.shape[-1]
    y = F.conv2d(h.unsqueeze(1).permute(0, 3, 1, 2), w.unsqueeze(2).to(memory_format=torch.channels_last),
                 padding=(0, d * (k // 2)), dilation=(1, d))
    return y.permute(0, 2, 3, 1).reshape(h.shape[0], h.shape[1], w.shape[0])
def bench(fn, n=20):
    for _ in range(3): fn()
    torch.cuda.synchronize(); s, e = torch.cuda.Event(True), torch.cuda.Event(True); s.record()
    for _ in range(n): fn()
    e.record(); torch.cuda.synchronize(); return s.elapsed_time(e) / n
print(torch.cuda.get_device_name())
specs = [(5, 1), (5, 2), (5, 4), (5, 8), (5, 16), (9, 1), (9, 4)]
print(f"{'L':>5s} {'C':>4s} " + " ".join(f"{f'k{k} d{d}':>9s}" for k, d in specs) + f" {'norm+k5 block':>14s}")
for L, C in [(1024, 48), (512, 64), (256, 80), (128, 96)]:
    x = torch.randn(N, L, C, device="cuda", dtype=torch.bfloat16, requires_grad=True)
    dy = torch.randn_like(x)
    row = []
    for k, d in specs:
        w = (torch.randn(C, C, k, device="cuda") * 0.02).to(torch.bfloat16).requires_grad_()
        row.append(f"{bench(lambda: torch.autograd.grad(conv(x, w, d), [x, w], dy)):9.2f}")
    g = torch.ones(C, device="cuda", requires_grad=True); b = torch.zeros(C, device="cuda", requires_grad=True)
    w = (torch.randn(C, C, 5, device="cuda") * 0.02).to(torch.bfloat16).requires_grad_()
    blk = bench(lambda: torch.autograd.grad(conv(read_norm_film_gelu(x, g, b, None, 32), w, 1), [x, w, g, b], dy))
    print(f"{L:5d} {C:4d} " + " ".join(row) + f" {blk:14.2f}", flush=True)
    del x, dy; torch.cuda.empty_cache()
