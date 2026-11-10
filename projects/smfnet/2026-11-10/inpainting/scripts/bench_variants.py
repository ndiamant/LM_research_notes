"""Training-step time (get_loss, backward, clip, AdamW; compiled, bf16) and peak
memory of U-Net variants. Usage: bench_variants.py <variant> [smf_length regions reads]
(default 1024 bp x 256 regions x 32 reads per region)."""
import sys, time
from types import SimpleNamespace
import torch
from smf_net.cnn import CNNConfig
from smf_net.data import SMFItem
from smf_net.model import SMFNet, SMFNetConfig

FAST = dict(dna_encoder_type="alphagenome", smf_body_type="unet", unet_attention_layers=2,
            unet_cross_attention=True, unet_channels_last=True, output_conv_matmul=True)
READ = dict(unet_norm="read", unet_fused_norm=True)
POS = dict(unet_norm="position")
VARIANTS = {
    "k9": dict(**READ, unet_kernel_size=9),
    "k9_linear": dict(**READ, unet_kernel_size=9, unet_upsample="linear"),
    "k9_dilated": dict(**READ, unet_kernel_size=9, unet_down_dilations=[4, 4, 2, 1, 1, 1, 1]),
    "k9_local_attn": dict(**READ, unet_kernel_size=9, unet_local_attention_levels=[2]),
    "k9_bottleneck16": dict(**READ, unet_kernel_size=9, unet_channels=[48, 64, 80, 96, 112]),
    "k9_read_masked": dict(unet_norm="read_masked", unet_kernel_size=9),
    "k9_window": dict(unet_norm="window", unet_kernel_size=9),
    "base": dict(**READ),
    "base_pos": dict(**POS),
    "p2": dict(**READ, unet_patch_size=2, unet_channels=[64, 80, 96, 112, 128, 144]),
    "p4_narrow": dict(**READ, unet_patch_size=4, unet_channels=[80, 96, 112, 128, 144]),
    "p4_mid": dict(**READ, unet_patch_size=4, unet_channels=[96, 128, 160, 192, 224]),
    "p4_wide": dict(**READ, unet_patch_size=4, unet_channels=[128, 160, 192, 224, 256]),
    "p4_wide_pos": dict(**POS, unet_patch_size=4, unet_channels=[128, 160, 192, 224, 256]),
    "p4_wide_h2": dict(**READ, unet_patch_size=4, unet_channels=[128, 160, 192, 224, 256],
                       unet_head_channels=32, unet_head_blocks=2),
    "base_k3": dict(**READ, unet_kernel_size=3),
    "base_k7": dict(**READ, unet_kernel_size=7),
    "base_k9": dict(**READ, unet_kernel_size=9),
    "base_k15": dict(**READ, unet_kernel_size=15),
    "base_e1": dict(**READ, unet_extra_blocks=1),
    "base_e2": dict(**READ, unet_extra_blocks=2),
    "base_k9_e1": dict(**READ, unet_kernel_size=9, unet_extra_blocks=1),
    "p4_wide_h2_k9": dict(**READ, unet_patch_size=4, unet_channels=[128, 160, 192, 224, 256],
                          unet_head_channels=32, unet_head_blocks=2, unet_kernel_size=9),
    "p4_wide_h2_e1": dict(**READ, unet_patch_size=4, unet_channels=[128, 160, 192, 224, 256],
                          unet_head_channels=32, unet_head_blocks=2, unet_extra_blocks=1),
    "base_k9_d1_u5": dict(**READ, unet_kernel_size=9, unet_down_blocks=1, unet_up_kernel_size=5),
    "base_k9_d1": dict(**READ, unet_kernel_size=9, unet_down_blocks=1),
    "base_k13_d1_u5": dict(**READ, unet_kernel_size=13, unet_down_blocks=1, unet_up_kernel_size=5),
    "p4_mid_pos": dict(**POS, unet_patch_size=4, unet_channels=[96, 128, 160, 192, 224]),
}
name = sys.argv[1]
smf_len = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
B = int(sys.argv[3]) if len(sys.argv) > 3 else 256
R = int(sys.argv[4]) if len(sys.argv) > 4 else 32
L = 2 * smf_len
torch.backends.cudnn.benchmark = True
cfg = CNNConfig(**FAST, **VARIANTS[name])
model = SMFNet(SMFNetConfig(cnn_config=cfg, input_length=L)).cuda()
model.smart_log = lambda *a, **k: None
model._trainer = SimpleNamespace(max_epochs=400, estimated_stepping_batches=182000)
opt = model.configure_optimizers()["optimizer"]
model.model.compile()
body = sum(p.numel() for n, p in model.model.named_parameters() if not n.startswith("dna_encoder"))
g = torch.Generator().manual_seed(0)
batch = SMFItem(smf=torch.randint(-1, 2, (B, R, L), generator=g).cuda(),
                dna=torch.nn.functional.one_hot(torch.randint(4, (B, L), generator=g), 4).transpose(1, 2).float().cuda())
def step():
    with torch.autocast("cuda", dtype=torch.bfloat16):
        loss = model.get_loss(batch, "train")
    loss.backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
    opt.step(); opt.zero_grad(set_to_none=True)
for _ in range(5): step()
torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
t = time.time()
for _ in range(30): step()
torch.cuda.synchronize()
ms = (time.time() - t) / 30 * 1e3
print(f"RESULT {name:12s} smf {smf_len:5d} bp x {B:3d} regions x {R:3d} reads: {ms:6.1f} ms/step, "
      f"{B * R / ms:5.1f}k reads/s, "
      f"peak {torch.cuda.max_memory_allocated() / 2**30:5.1f} GiB, non-DNA-encoder params {body / 1e6:.2f}M", flush=True)
