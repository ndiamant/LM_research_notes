"""Full-size training-step time and peak memory: baseline CNN vs Cheri blocks."""
import sys
import time

import torch

from smf_net.model import SMFNet, SMFNetConfig
from smf_net.cnn import CNNConfig

torch.manual_seed(0)
R, L = 32, 2048
variants = {
    "cnn baseline (AdamW), B=256": (256, dict(), "adamw"),
    "cheri + muon, B=256": (256, dict(block_type="cheri"), "muon"),
    "cheri + muon + checkpoint, B=256": (256, dict(block_type="cheri", cheri_checkpoint=True), "muon"),
    "cheri + muon, B=128": (128, dict(block_type="cheri"), "muon"),
}
only = sys.argv[1:] or list(variants)
for name in only:
    B, cnn_kw, opt = variants[name]
    model = optimizer = None
    try:
        model = SMFNet(SMFNetConfig(cnn_config=CNNConfig(**cnn_kw), optimizer=opt)).cuda()
        from types import SimpleNamespace
        model._trainer = SimpleNamespace(max_epochs=400, estimated_stepping_batches=182000)
        optimizer = model.configure_optimizers()["optimizer"]
        n_params = sum(p.numel() for p in model.parameters()) / 1e6
        model.model.compile()
        dna = torch.nn.functional.one_hot(torch.randint(4, (B, L)), 4).transpose(1, 2).float().cuda()
        smf = torch.randint(0, 3, (B, R, 3, L // 2)).float().cuda()
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
        per_region = dt / B * 256
        print(f"{name}: params {n_params:.1f}M, step {dt*1000:.0f} ms ({per_region*1000:.0f} ms per 256 regions), "
              f"peak {torch.cuda.max_memory_allocated()/2**30:.1f} GiB, 400 epochs ~ {per_region*455*400/3600:.1f} h", flush=True)
    except torch.OutOfMemoryError:
        print(f"{name}: out of memory", flush=True)
    del model, optimizer
    torch.cuda.empty_cache()
