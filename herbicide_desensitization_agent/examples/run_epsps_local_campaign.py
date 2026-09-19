"""Run the available EPSPS campaign on a prepared Linux Docker/GPU host."""
from __future__ import annotations

import argparse
import importlib.metadata
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

from .epsps_local_campaign import collect, prepare, write_json

IMAGE = "rbgcsail/diffdock@sha256:1b7bb3adb332fdc9648a0ec53dec2f790cfbb816d7478d2bd93f1cdea3b269f0"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "workflow", "contacts", "cached", "output", "diffdock-cache", "boltz-cache"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--git-commit", required=True)
    args = parser.parse_args()
    if importlib.metadata.version("boltz") != "2.2.1":
        raise ValueError("This protocol requires boltz==2.2.1")
    root = args.output.resolve()
    counts = prepare(args.inputs, args.workflow, args.contacts, args.cached, root)
    print(counts, flush=True)
    (root / "logs").mkdir()
    script = Path(__file__).with_name("run_diffdock_smoke.py")
    import shutil
    shutil.copy2(script, root / "run_diffdock.py")
    state = {"git_commit": args.git_commit, "started_utc": datetime.now(timezone.utc).isoformat(),
             "status": "RUNNING", "counts": counts, "stages": [],
             "blocked_stages": {
                 "five_other_targets": "Validated live complex inputs are not configured.",
                 "evidence_agents_and_llm_judge": "Live integration/access is not configured; no mock reviews used.",
                 "kd_equivalence_and_fold_ddg": "Matched Kd intervals and folding ddG unavailable; Boltz pIC50 is not Kd.",
                 "experimental_learning": "No real assays supplied.",
             }}

    def save():
        write_json(root / "execution_manifest.json", state)

    def start(name, command, env=None):
        row = {"stage": name, "command": command, "status": "RUNNING",
               "started_utc": datetime.now(timezone.utc).isoformat()}
        state["stages"].append(row)
        save()
        stream = (root / "logs" / f"{name}.log").open("w")
        try:
            process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
        except BaseException:
            stream.close()
            raise
        return process, stream, row

    def finish(job):
        process, stream, row = job
        code = process.wait()
        stream.close()
        row.update(status="COMPLETED" if code == 0 else "FAILED", exit_code=code,
                   finished_utc=datetime.now(timezone.utc).isoformat())
        save()
        print(f"{row['stage']}: {row['status']}", flush=True)
        if code:
            raise RuntimeError(f"{row['stage']} failed; see logs/{row['stage']}.log")

    jobs, containers = [], []
    try:
        for condition in ("native", "herbicide"):
            for replicate in (1, 2):
                shard = f"{condition}-{replicate}"
                container = f"herbicide-campaign-{shard}"
                command = ["docker", "run", "--rm", "--name", container, "--network", "none", "--shm-size=4g",
                           "-e", "PYTHONPATH=/home/appuser/DiffDock", "-e", "OMP_NUM_THREADS=2",
                           "-e", "DIFFDOCK_THREADS=2", "-e", "DIFFDOCK_DEVICE=cpu",
                           "-e", f"DIFFDOCK_SEED={42 + replicate}",
                           "-e", f"DIFFDOCK_CSV=/validation/inputs/diffdock-{shard}.csv",
                           "-e", f"DIFFDOCK_OUTPUT=/validation/diffdock/{shard}",
                           "-v", f"{root}:/validation",
                           "-v", f"{args.diffdock_cache.resolve()}/workdir:/home/appuser/DiffDock/workdir:ro",
                           "-v", f"{args.diffdock_cache.resolve()}/torch:/home/appuser/.cache/torch:ro",
                           IMAGE, "micromamba", "run", "-n", "diffdock", "python", "-u", "/validation/run_diffdock.py"]
                jobs.append(start(f"diffdock-{shard}", command))
                containers.append(container)
        env = {**os.environ, "OMP_NUM_THREADS": "4", "MKL_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "4"}
        for replicate in (1, 2):
            command = [str(Path(sys.executable).with_name("boltz")), "predict",
                       str(root / "inputs" / f"replicate-{replicate}"), "--out_dir",
                       str(root / "boltz" / f"replicate-{replicate}"), "--cache", str(args.boltz_cache.resolve()),
                       "--model", "boltz2", "--devices", "1", "--accelerator", "gpu",
                       "--recycling_steps", "3", "--sampling_steps", "100", "--diffusion_samples", "1",
                       "--diffusion_samples_affinity", "3", "--sampling_steps_affinity", "200",
                       "--seed", str(42 + replicate), "--num_workers", "0", "--no_kernels"]
            job = start(f"boltz-replicate-{replicate}", command, env)
            jobs.append(job)
            finish(job)
        for job in jobs:
            if job[2]["status"] == "RUNNING":
                finish(job)
        state["validation"] = collect(root)
        state["status"] = "COMPUTATION_COMPLETE_ANALYSIS_PENDING"
        print(state["validation"], flush=True)
    except BaseException:
        state["status"] = "FAILED"
        raise
    finally:
        for process, stream, row in jobs:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                row.update(status="INTERRUPTED", exit_code=process.returncode)
            stream.close()
        if state["status"] == "FAILED":
            for container in containers:
                subprocess.run(["docker", "stop", "--time", "5", container], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        state["finished_utc"] = datetime.now(timezone.utc).isoformat()
        save()


if __name__ == "__main__":
    main()
