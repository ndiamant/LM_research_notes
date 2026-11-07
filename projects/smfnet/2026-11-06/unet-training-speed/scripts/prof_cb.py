"""Lightning callback: during training steps 60-69, warn on every host-device
sync, time each step on the host, and record a CPU+CUDA trace."""
import time, traceback, warnings
import torch
from pytorch_lightning import Callback
from torch.profiler import ProfilerActivity, profile

def _showwarning(message, category, filename, lineno, file=None, line=None):
    stack = "".join(traceback.format_stack(limit=12)[:-2])
    print(f"SYNC: {message}\n{stack}", flush=True)

class ProfileSteps(Callback):
    t = None
    def on_train_batch_start(self, trainer, pl_module, batch, batch_idx):
        now = time.perf_counter()
        if self.t is not None and 50 <= batch_idx <= 75:
            print(f"STEP {batch_idx}: {1e3 * (now - self.t):.0f} ms", flush=True)
        self.t = now
        if batch_idx == 60:
            warnings.showwarning = _showwarning
            torch.cuda.set_sync_debug_mode("warn")
        if batch_idx == 64:
            torch.cuda.set_sync_debug_mode(0)
            torch.cuda.synchronize()
            self.prof = profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA]); self.prof.__enter__()
    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        if batch_idx == 69:
            torch.cuda.synchronize(); self.prof.__exit__(None, None, None)
            self.prof.export_chrome_trace("/scratch/users/diamant/smfnet_fast/lightning_trace.json")
            print("PROFILE: " + self.prof.key_averages().table(sort_by="self_cpu_time_total", row_limit=25).replace("\n", "\nPROFILE: "))
