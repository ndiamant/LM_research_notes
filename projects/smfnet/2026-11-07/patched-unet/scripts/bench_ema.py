"""Time one EMA update of the whole SMFNet (as Lightning's EMAWeightAveraging
does: AveragedModel over all parameters and buffers): per-tensor avg_fn (the
callback's default) vs foreach multi_avg_fn."""
import time
import torch
from torch.optim.swa_utils import AveragedModel, get_ema_avg_fn, get_ema_multi_avg_fn
from smf_net.cnn import CNNConfig
from smf_net.model import SMFNet, SMFNetConfig
cfg = CNNConfig(dna_encoder_type="alphagenome", smf_body_type="unet", unet_norm="read", unet_fused_norm=True,
                unet_attention_layers=2, unet_cross_attention=True, unet_channels_last=True, output_conv_matmul=True,
                unet_patch_size=4, unet_channels=[128, 160, 192, 224, 256], unet_head_channels=32, unet_head_blocks=2)
model = SMFNet(SMFNetConfig(cnn_config=cfg)).cuda()
n = sum(p.numel() for p in model.parameters()); t = len(list(model.parameters())) + len(list(model.buffers()))
print(f"{n / 1e6:.0f}M parameters, {t} tensors")
for name, kw in [("per-tensor avg_fn", dict(avg_fn=get_ema_avg_fn(0.999))),
                 ("foreach multi_avg_fn", dict(multi_avg_fn=get_ema_multi_avg_fn(0.999)))]:
    ema = AveragedModel(model, use_buffers=True, **kw)
    for _ in range(5): ema.update_parameters(model)
    torch.cuda.synchronize(); s = time.time()
    for _ in range(20): ema.update_parameters(model)
    torch.cuda.synchronize(); print(f"{name}: {(time.time() - s) / 20 * 1e3:.1f} ms per update")
