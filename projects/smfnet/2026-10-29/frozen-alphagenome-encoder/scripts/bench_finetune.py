import sys, time, torch
from smf_net.cnn import CNN, CNNConfig
from smf_net.alphagenome import load_alphagenome_encoder_weights

W = sys.argv[1]
torch.manual_seed(0)
B, R, L = 256, 32, 2048
dna = torch.nn.functional.one_hot(torch.randint(4, (B, L)), 4).transpose(1, 2).float().cuda()
smf = torch.randn(B, R, 3, L // 2).cuda()
variants = {
    "cnn baseline": dict(dna_encoder_type="cnn"),
    "AG frozen": dict(dna_encoder_type="alphagenome"),
    "AG train last 3 blocks": dict(dna_encoder_type="alphagenome", alphagenome_trainable_blocks=3),
    "AG train all": dict(dna_encoder_type="alphagenome", alphagenome_frozen=False),
}
for name, kw in variants.items():
    model = CNN(CNNConfig(**kw)).cuda()
    if kw["dna_encoder_type"] == "alphagenome":
        load_alphagenome_encoder_weights(model.dna_encoder.alphagenome, W)
    nt = sum(p.numel() for p in model.parameters() if p.requires_grad) / 1e6
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad], lr=1e-3)
    model.compile()
    torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
    def step():
        with torch.autocast("cuda", dtype=torch.bfloat16):
            loss = model(smf, dna).float().square().mean()
        loss.backward(); opt.step(); opt.zero_grad(set_to_none=True)
    for _ in range(5): step()
    torch.cuda.synchronize(); t = time.time()
    for _ in range(20): step()
    torch.cuda.synchronize(); dt = (time.time() - t) / 20
    print(f"{name}: trainable {nt:.1f}M, step {dt*1000:.0f} ms, peak {torch.cuda.max_memory_allocated()/2**30:.1f} GiB, "
          f"400 epochs x 455 steps ~ {dt*455*400/3600:.1f} h", flush=True)
    del model, opt; torch.cuda.empty_cache()
