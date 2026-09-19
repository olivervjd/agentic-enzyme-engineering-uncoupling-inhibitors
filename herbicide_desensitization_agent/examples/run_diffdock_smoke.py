"""Run inside the official DiffDock-L container with /validation mounted."""

import io
import csv
import json
import os
import random
from pathlib import Path

# This pinned image's old CUDA/PyTorch build does not support H200 (sm_90).
# Default to CPU; opt in to GPU only with a compatible runtime.
if os.environ.get("DIFFDOCK_DEVICE", "cpu") == "cpu":
    os.environ["CUDA_VISIBLE_DEVICES"] = ""

import numpy as np
import torch
import yaml

import inference


def check_graph(graph, stage):
    for entity in ("ligand", "receptor"):
        for field in ("pos", "x"):
            values = graph[entity][field]
            if values.numel() == 0 or not torch.isfinite(values).all():
                invalid = (~torch.isfinite(values)).nonzero().cpu().tolist()[:5]
                raise ValueError(f"{stage}: empty or non-finite {entity}.{field}, "
                                 f"shape={tuple(values.shape)}, first invalid indices={invalid}")


original_randomize = inference.randomize_position


def checked_randomize(graphs, *args, **kwargs):
    for graph in graphs:
        check_graph(graph, "before randomization")
    result = original_randomize(graphs, *args, **kwargs)
    for graph in graphs:
        check_graph(graph, "after randomization")
    return result


original_get_model = inference.get_model


def checked_get_model(*args, **kwargs):
    model = original_get_model(*args, **kwargs)

    def check_output(module, inputs, outputs):
        for index, value in enumerate(outputs if isinstance(outputs, tuple) else (outputs,)):
            if isinstance(value, torch.Tensor) and not torch.isfinite(value).all():
                raise ValueError(f"Non-finite DiffDock model output {index}")

    model.register_forward_hook(check_output)
    return model


inference.randomize_position = checked_randomize
inference.get_model = checked_get_model
seed = int(os.environ.get("DIFFDOCK_SEED", "42"))
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.set_num_threads(int(os.environ.get("DIFFDOCK_THREADS", "4")))
csv_path = Path(os.environ.get("DIFFDOCK_CSV", "/validation/inputs/diffdock.csv"))
out = Path(os.environ.get("DIFFDOCK_OUTPUT", "/validation/diffdock"))
with csv_path.open() as stream:
    names = [row["complex_name"] for row in csv.DictReader(stream)]
if not names or len(set(names)) != len(names) or any(Path(name).name != name for name in names):
    raise ValueError("Expected unique, simple complex names")
args = inference.get_parser().parse_args([
    "--protein_ligand_csv", str(csv_path),
    "--out_dir", str(out), "--batch_size", "4", "--loglevel", "INFO",
])
config = yaml.safe_load(args.config)
args.config.close()
# Upstream YAML overrides CLI flags; set sample count in the configuration itself.
config["samples_per_complex"] = 4
args.config = io.StringIO(yaml.safe_dump(config))
for name in names:
    if list((out / name).glob("rank*.sdf")):
        raise FileExistsError("Use a fresh output directory; refusing to validate stale docking poses")
out.mkdir(parents=True, exist_ok=True)
# Preserve the original smoke-test manifest paths while isolating campaign shards.
metadata_dir = out.parent if out == Path("/validation/diffdock") else out
(metadata_dir / "diffdock_config.json").write_text(json.dumps(config, indent=2))
(metadata_dir / "diffdock_runtime.json").write_text(json.dumps({
    "torch_version": torch.__version__, "cuda_version": torch.version.cuda,
    "device": "cuda" if torch.cuda.is_available() else "cpu", "seed": seed,
}, indent=2))
inference.main(args)
for name in names:
    poses = list((out / name).glob("rank*_confidence*.sdf"))
    if len(poses) != 4:
        raise RuntimeError(f"DiffDock did not produce four ranked poses for {name}")
print(f"DIFFDOCK_EXECUTION_PASS: {4 * len(names)} ranked poses across {len(names)} complexes", flush=True)
