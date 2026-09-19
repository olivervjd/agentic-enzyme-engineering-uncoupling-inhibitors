"""Run inside the official DiffDock-L container with /validation mounted."""

import io
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
random.seed(42)
np.random.seed(42)
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
torch.set_num_threads(4)
args = inference.get_parser().parse_args([
    "--protein_ligand_csv", "/validation/inputs/diffdock.csv",
    "--out_dir", "/validation/diffdock", "--batch_size", "4", "--loglevel", "INFO",
])
config = yaml.safe_load(args.config)
args.config.close()
# Upstream YAML overrides CLI flags; set sample count in the configuration itself.
config["samples_per_complex"] = 4
args.config = io.StringIO(yaml.safe_dump(config))
for ligand in ("pep", "glyphosate"):
    if list((Path("/validation/diffdock") / f"epsps-{ligand}").glob("rank*.sdf")):
        raise FileExistsError("Use a fresh output directory; refusing to validate stale docking poses")
Path("/validation/diffdock_config.json").write_text(json.dumps(config, indent=2))
Path("/validation/diffdock_runtime.json").write_text(json.dumps({
    "torch_version": torch.__version__, "cuda_version": torch.version.cuda,
    "device": "cuda" if torch.cuda.is_available() else "cpu", "seed": 42,
}, indent=2))
inference.main(args)
for ligand in ("pep", "glyphosate"):
    poses = list((Path("/validation/diffdock") / f"epsps-{ligand}").glob("rank*_confidence*.sdf"))
    if len(poses) != 4:
        raise RuntimeError(f"DiffDock did not produce four ranked poses for {ligand}")
print("DIFFDOCK_EXECUTION_PASS: 8 ranked poses across 2 WT ligand complexes", flush=True)
