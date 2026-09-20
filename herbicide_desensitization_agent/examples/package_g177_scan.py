"""Package an explicit allowlist of scientific scan artifacts; never credentials."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument("root", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
root = args.root.resolve()
results = json.loads((root / "results.json").read_text())
validation = json.loads((root / "validation.json").read_text())
manifest = json.loads((root / "manifest.json").read_text())
if not results["complete"] or validation["status"] != "PASSED" or validation["validated_jobs"] != 120:
    raise ValueError("Scientific outputs are not all validated")
names = {"manifest.json", "status.json", "execution.json", "results.json", "validation.json", "variants.fasta", "binding_table.md"}
for record in manifest["records"]:
    names.update((record["request"], record["msa"]))
for record in results["replicates"]:
    for artifact in record["artifacts"]:
        path = (root / artifact["path"]).resolve()
        path.relative_to(root)
        if hashlib.sha256(path.read_bytes()).hexdigest() != artifact["sha256"]:
            raise ValueError("Artifact hash changed")
        names.add(artifact["path"])
if args.output.exists():
    raise FileExistsError("Refusing to overwrite archive")
with tarfile.open(args.output, "w:gz") as archive:
    for name in sorted(names):
        target = (root / name).resolve()
        target.relative_to(root)
        if not target.is_file():
            raise ValueError("Missing allowlisted artifact")
        archive.add(target, arcname="g177-scan/" + name, recursive=False)
print(json.dumps({"archive": str(args.output), "file_count": len(names), "bytes": args.output.stat().st_size,
    "sha256": hashlib.sha256(args.output.read_bytes()).hexdigest()}))
