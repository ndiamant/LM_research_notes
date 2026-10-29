import sys, time, torch
from smf_net.cnn import CNN, CNNConfig
from smf_net.alphagenome import load_alphagenome_encoder_weights

W = sys.argv[1]
torch.manual_seed(0)
B, R, L = 256, 32, 2048
dna = torch.nn.functional.one_hot(torch.randint(4, (B, L)), 4).transpose(1, 2).float().cuda()
smf = torch.randn(B, R, 3, L // 2).cuda()
for kind in ["cnn", "alphagenome"]:
    model = CNN(CNNConfig(dna_encoder_type=kind)).cuda()
    if kind == "alphagenome":
        load_alphagenome_encoder_weights(model.dna_encoder.alphagenome, W)
    n = sum(p.numel() for p in model.parameters()) / 1e6
    nt = sum(p.numel() for p in model.parameters() if p.requires_grad) / 1e6
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=1e-3)
    model.compile()
    torch.cuda.reset_peak_memory_stats()
    def step():
        with torch.autocast("cuda", dtype=torch.bfloat16):
            loss = model(smf, dna).float().square().mean()
        loss.backward(); opt.step(); opt.zero_grad(set_to_none=True)
        return loss
    for _ in range(5): step()
    torch.cuda.synchronize(); t = time.time()
    for _ in range(20): loss = step()
    torch.cuda.synchronize()
    dt = (time.time() - t) / 20
    with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
        emb = model.dna_encoder(dna[:8]).float()
        feats = model.dna_encoder.alphagenome(dna[:8]) if kind == "alphagenome" else []
    print(f"{kind}: params {n:.1f}M trainable {nt:.1f}M step {dt*1000:.0f} ms "
          f"peak {torch.cuda.max_memory_allocated()/2**30:.1f} GiB loss {loss.item():.3f} "
          f"embed rms {emb.square().mean().sqrt():.2f} finite {torch.isfinite(emb).all().item()}")
    for f in feats: print("   AG feature rms", f.shape[1], f"{f.float().square().mean().sqrt():.2f}")
    del model, opt; torch.cuda.empty_cache()
