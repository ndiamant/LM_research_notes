"""Full-size training-step time and peak memory: frozen-AG encoder with the
dilated SMF body vs U-Net SMF bodies."""
import time
from types import SimpleNamespace

import torch

from smf_net.cnn import CNNConfig
from smf_net.model import SMFNet, SMFNetConfig

torch.manual_seed(0)
B, R, L = 256, 32, 2048
variants = {
    "AG frozen + dilated body (current)": dict(),
    "AG frozen + U-Net 48..144 (same params)": dict(smf_body_type="unet"),
    "AG frozen + U-Net 64,128..384 (wide)": dict(smf_body_type="unet", unet_channels=[64, 128, 160, 224, 288, 352, 384]),
}
dna = torch.nn.functional.one_hot(torch.randint(4, (B, L)), 4).transpose(1, 2).float().cuda()
smf = torch.randint(0, 3, (B, R, 3, L // 2)).float().cuda()
for name, kw in variants.items():
    model = optimizer = None
    try:
        model = SMFNet(SMFNetConfig(cnn_config=CNNConfig(dna_encoder_type="alphagenome", **kw))).cuda()
        model._trainer = SimpleNamespace(max_epochs=400, estimated_stepping_batches=182000)
        optimizer = model.configure_optimizers()["optimizer"]
        body = sum(p.numel() for p in model.model.body.parameters()) / 1e6
        trainable = sum(p.numel() for p in model.parameters() if p.requires_grad) / 1e6
        model.model.compile()
        torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
        def step():
            with torch.autocast("cuda", dtype=torch.bfloat16):
                loss = model.model(smf, dna).float().square().mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step(); optimizer.zero_grad(set_to_none=True)
        for _ in range(5): step()
        torch.cuda.synchronize(); t = time.time()
        for _ in range(15): step()
        torch.cuda.synchronize(); dt = (time.time() - t) / 15
        print(f"{name}: SMF body {body:.2f}M params, trainable {trainable:.1f}M, step {dt*1000:.0f} ms, "
              f"peak {torch.cuda.max_memory_allocated()/2**30:.1f} GiB, 400 epochs ~ {dt*455*400/3600:.1f} h", flush=True)
    except torch.OutOfMemoryError:
        print(f"{name}: out of memory", flush=True)
    del model, optimizer
    torch.cuda.empty_cache()
