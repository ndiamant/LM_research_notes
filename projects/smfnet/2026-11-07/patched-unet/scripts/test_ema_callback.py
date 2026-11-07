"""ForeachEMAWeightAveraging: same average as Lightning's EMAWeightAveraging,
and it loads its state from a run saved with Lightning's callback."""
import sys, torch
sys.path.insert(0, "/home/users/diamant/repos/SMFNet_patch/scripts")
from pytorch_lightning.callbacks import EMAWeightAveraging
from torch.optim.swa_utils import AveragedModel
from train import ForeachEMAWeightAveraging
torch.manual_seed(0)
m = torch.nn.Sequential(torch.nn.Linear(4, 4), torch.nn.BatchNorm1d(4)).cuda()
a, b = EMAWeightAveraging(decay=0.9), ForeachEMAWeightAveraging(decay=0.9)
assert a.state_key == b.state_key
for cb in (a, b):
    cb._average_model = AveragedModel(m, use_buffers=True, **cb._kwargs)
for _ in range(5):
    with torch.no_grad():
        for p in m.parameters(): p.add_(torch.randn_like(p))
    m(torch.randn(8, 4, device="cuda"))
    a._average_model.update_parameters(m); b._average_model.update_parameters(m)
for x, y in zip(a._average_model.state_dict().values(), b._average_model.state_dict().values()):
    torch.testing.assert_close(x, y)
ck = torch.load("/scratch/users/diamant/smf_models/ONT_mESC_unet_attn_adamw_ema_p4_narrow/last.ckpt", map_location="cpu", weights_only=False)
print("checkpoint callback keys:", list(ck["callbacks"]))
print("ok")
