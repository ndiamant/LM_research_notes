"""How far the SMF predictor sees: gradient magnitude by distance from one output."""

import torch

from smf_net.cnn import CNN, CNNConfig

UNET = dict(
    dna_encoder_type="alphagenome", smf_body_type="unet", unet_kernel_size=9,
    unet_down_dilations=[4, 4, 2, 1, 1, 1, 1], unet_channels_last=True, output_conv_matmul=True,
)


def profile(norm, attention_layers, smf_len=8192, flank=1024, seeds=3):
    extra = {"unet_norm": norm, "unet_attention_layers": attention_layers,
             "unet_cross_attention": attention_layers > 0, "dna_flank": flank,
             "unet_fused_norm": norm == "read"}
    total = torch.zeros(smf_len, device="cuda")
    for seed in range(seeds):
        torch.manual_seed(seed)
        model = CNN(CNNConfig(**(UNET | extra))).cuda().float()
        with torch.no_grad():
            for p in model.parameters():
                p.normal_(0, 0.05 if p.ndim > 1 else 0.5)
        smf = torch.randn(1, 1, 3, smf_len, device="cuda", requires_grad=True)
        dna = torch.nn.functional.one_hot(torch.randint(4, (1, smf_len + 2 * flank)), 4)
        model(smf, dna.transpose(1, 2).float().cuda())[0, 0, 0, smf_len // 2].backward()
        total += smf.grad[0, 0].abs().sum(0)
    g = total / total.max()
    d = (torch.arange(smf_len, device="cuda") - smf_len // 2).abs()
    out = {}
    for thr in (1e-2, 1e-3, 1e-4):
        out[thr] = int(d[g > thr].max())
    support = int(d[g > 0].max())
    return support, out


for norm, layers in (("position", 0), ("position", 2), ("read", 0), ("read", 2)):
    support, reach = profile(norm, layers)
    print(f"norm={norm:8s} attention={layers}: nonzero to +/-{support} bp; gradient above "
          + ", ".join(f"{t:g} of peak to +/-{r} bp" for t, r in reach.items()))
