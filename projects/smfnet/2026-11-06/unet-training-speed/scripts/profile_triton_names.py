"""Profile training steps of the attention U-Net (frozen AG encoder), full
batch (256 regions x 32 reads), bf16 autocast. mode=compiled (as in training):
CUDA time by kernel and kernel category. mode=eager: CUDA time by aten op,
plus by module for attribution."""
import sys
from collections import defaultdict
from types import SimpleNamespace

import torch
from torch.profiler import ProfilerActivity, profile, record_function

from smf_net.cnn import CNNConfig
from smf_net.model import SMFNet, SMFNetConfig

mode = sys.argv[1]
import torch._inductor.config as icfg
icfg.triton.unique_kernel_names = True
icfg.triton.descriptive_names = 'original_aten'
torch.backends.cudnn.benchmark = True
B, R, L = 256, 32, 2048
cfg = CNNConfig(dna_encoder_type="alphagenome", smf_body_type="unet", unet_norm="read",
                unet_attention_layers=2, unet_cross_attention=True, unet_fused_norm="--fused" in sys.argv, unet_channels_last="--cl" in sys.argv, output_conv_matmul="--cl" in sys.argv)
model = SMFNet(SMFNetConfig(cnn_config=cfg)).cuda()
model._trainer = SimpleNamespace(max_epochs=400, estimated_stepping_batches=182000)
opt = model.configure_optimizers()["optimizer"]
if mode == "compiled":
    model.model.compile()
dna = torch.nn.functional.one_hot(torch.randint(4, (B, L)), 4).transpose(1, 2).float().cuda()
smf = torch.nn.functional.one_hot(torch.randint(0, 3, (B, R, L // 2)), 3).permute(0, 1, 3, 2).float().cuda()

def step():
    with torch.autocast("cuda", dtype=torch.bfloat16):
        with record_function("forward"):
            loss = model.model(smf, dna).float().square().mean()
    with record_function("backward"):
        loss.backward()
    with record_function("optimizer"):
        opt.step(); opt.zero_grad(set_to_none=True)

for _ in range(5): step()
torch.cuda.synchronize()
hooks = []
if mode == "eager":  # attribute CUDA time to the model's main parts
    parts = {"dna_encoder.alphagenome": model.model.dna_encoder.alphagenome,
             "dna_encoder.projections": model.model.dna_encoder.projections,
             "body.down": model.model.body.down, "body.attention": model.model.body.attention,
             "body.up": model.model.body.up}
    for name, mod in parts.items():
        ctx = {}
        hooks.append(mod.register_forward_pre_hook(lambda m, a, n=name, c=ctx: c.__setitem__("r", record_function("part:" + n).__enter__())))
        hooks.append(mod.register_forward_hook(lambda m, a, o, c=ctx: c["r"].__exit__(None, None, None)))
with profile(activities=[ProfilerActivity.CUDA, ProfilerActivity.CPU]) as prof:
    for _ in range(3): step()
    torch.cuda.synchronize()
events = prof.key_averages()
ANNOTATIONS = ("forward", "backward", "optimizer")
is_kernel = lambda e: e.device_type == torch.autograd.DeviceType.CUDA and e.key not in ANNOTATIONS and not e.key.startswith("part:")
kernel_time = sum(e.self_device_time_total for e in events if is_kernel(e))
print(f"{mode}: GPU kernel time per step {kernel_time / 3 / 1e3:.0f} ms")

def category(name):
    n = name.lower()
    if "nchwtonhwc" in n or "nhwctonchw" in n:
        return "cuDNN layout transposes (NCHW<->NHWC)"
    if any(k in n for k in ("gemm", "cutlass", "sm90", "sm80", "xmma", "cudnn", "conv", "wgrad", "dgrad", "fprop", "implicit")):
        return "matmul/conv"
    if "flash" in n or "fmha" in n or "attention" in n:
        return "attention"
    if "_stats_kernel" in n or "_forward_kernel" in n or "_backward_kernel" in n or "_backward_stats" in n:
        return "fused norm kernels"
    if "triton" in n and ("red" in n or "per_" in n):
        return "triton reduction (norms, sums)"
    if "triton" in n:
        return "triton pointwise (fused elementwise)"
    if "elementwise" in n or "vectorized" in n or "unrolled" in n:
        return "aten elementwise"
    if "reduce" in n or "norm" in n:
        return "aten reduction/norm"
    if "nchwtonhwc" in n or "nhwctonchw" in n:
        return "cuDNN layout transposes (NCHW<->NHWC)"
    if "copy" in n or "memcpy" in n or "memset" in n or "transpose" in n:
        return "copies/transposes"
    return "other"

if mode == "compiled":
    cats = defaultdict(float)
    for e in events:
        if is_kernel(e):
            cats[category(e.key)] += e.self_device_time_total
    print("\nby kernel category (share of GPU kernel time):")
    for c, t in sorted(cats.items(), key=lambda kv: -kv[1]):
        print(f"  {c:40s} {t / kernel_time:6.1%}  {t / 3 / 1e3:6.1f} ms/step")
    print("\ntop compiler-generated kernels:")
    tri = sorted((e for e in events if is_kernel(e) and e.key.startswith("triton_")), key=lambda e: -e.self_device_time_total)
    for e in tri[:30]:
        print(f"  {e.self_device_time_total / 3 / 1e3:6.2f} ms  {e.count // 3:3d}x  {e.key[:120]}")
    print("\ntop kernels:")
    for e in sorted((e for e in events if is_kernel(e)),
                    key=lambda e: -e.self_device_time_total)[:25]:
        print(f"  {e.self_device_time_total / kernel_time:6.1%}  {e.count // 3:4d}x  {e.key[:110]}")
else:
    print("\nby aten op (self CUDA time):")
    rows = sorted((e for e in events if e.key.startswith("aten::") and e.self_device_time_total > 0),
                  key=lambda e: -e.self_device_time_total)[:22]
    for e in rows:
        print(f"  {e.self_device_time_total / kernel_time:6.1%}  {e.count // 3:5d}x  {e.key}")
    print("\nby part (forward only, CUDA time incl. children):")
    for e in events:
        if e.key.startswith("part:"):
            print(f"  {e.key:32s} {e.device_time_total / 3 / 1e3:7.1f} ms/step")
    for e in events:
        if e.key in ("forward", "backward", "optimizer"):
            print(f"  {e.key:32s} {e.device_time_total / 3 / 1e3:7.1f} ms/step")
