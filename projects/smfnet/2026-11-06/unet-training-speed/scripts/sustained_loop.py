"""Full training step (get_loss, backward, clip, AdamW) without Lightning, with
and without a per-step .item() (as per-step logging does)."""
import time
from types import SimpleNamespace
import torch
from smf_net.cnn import CNNConfig
from smf_net.data import SMFItem
from smf_net.model import SMFNet, SMFNetConfig
torch.backends.cudnn.benchmark = True
B, R, L = 256, 32, 2048
cfg = CNNConfig(dna_encoder_type="alphagenome", smf_body_type="unet", unet_norm="read", unet_attention_layers=2,
                unet_cross_attention=True, unet_fused_norm=True, unet_channels_last=True, output_conv_matmul=True)
model = SMFNet(SMFNetConfig(cnn_config=cfg)).cuda()
model.smart_log = lambda *a, **k: None
model._trainer = SimpleNamespace(max_epochs=400, estimated_stepping_batches=182000)
opt = model.configure_optimizers()["optimizer"]
model.model.compile()
g = torch.Generator().manual_seed(0)
batch = SMFItem(smf=torch.randint(-1, 2, (B, R, L), generator=g).cuda(),
                dna=torch.nn.functional.one_hot(torch.randint(4, (B, L), generator=g), 4).transpose(1, 2).float().cuda())
def step(sync):
    with torch.autocast("cuda", dtype=torch.bfloat16):
        loss = model.get_loss(batch, "train")
    loss.backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
    opt.step(); opt.zero_grad(set_to_none=True)
    if sync:
        loss.item()
import subprocess
for _ in range(5): step(False)
torch.cuda.synchronize()
for chunk in range(10):
    t = time.time()
    for _ in range(20): step(False)
    torch.cuda.synchronize()
    clock = subprocess.run(["nvidia-smi", "--query-gpu=clocks.sm,power.draw,temperature.gpu", "--format=csv,noheader"], capture_output=True, text=True).stdout.strip()
    print(f"  steps {chunk*20:3d}-{chunk*20+19}: {(time.time() - t) / 20 * 1000:.0f} ms per step   [{clock}]", flush=True)
