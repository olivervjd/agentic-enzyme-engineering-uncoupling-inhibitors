from __future__ import annotations

import shutil
import math
import subprocess
import tempfile
from pathlib import Path


class FoldseekStructuralMatcher:
    """Pairwise Foldseek/TM-align adapter for mutant-versus-reference structures."""

    def __init__(self, executable: str = "foldseek") -> None:
        resolved = shutil.which(executable)
        if not resolved:
            raise RuntimeError("Foldseek executable is not installed")
        self.executable = resolved

    def compare(self, mutant: str | Path, reference: str | Path) -> dict[str, float]:
        mutant, reference = Path(mutant).resolve(), Path(reference).resolve()
        if not mutant.is_file() or not reference.is_file():
            raise FileNotFoundError("Both mutant and reference structures are required")
        with tempfile.TemporaryDirectory(prefix="foldseek-") as directory:
            result = Path(directory) / "alignment.tsv"
            subprocess.run([
                self.executable, "easy-search", str(mutant), str(reference), str(result),
                str(Path(directory) / "tmp"), "--alignment-type", "1", "--format-output",
                "query,target,qtmscore,ttmscore,lddt,qcov,tcov", "-a", "1",
            ], check=True, capture_output=True, text=True)
            lines = [line for line in result.read_text(encoding="utf-8").splitlines() if line.strip()]
            if len(lines) != 1:
                raise ValueError("Expected exactly one Foldseek alignment; select one protein chain explicitly")
            line = lines[0].split("\t")
            if len(line) != 7:
                raise ValueError("Unexpected Foldseek output schema")
        metrics = {
            "query_tm_score": float(line[2]), "target_tm_score": float(line[3]),
            "alignment_lddt": float(line[4]), "alignment_coverage": min(float(line[5]), float(line[6])),
        }
        if any(not math.isfinite(value) or not 0 <= value <= 1 for value in metrics.values()):
            raise ValueError("Invalid Foldseek metric")
        return metrics
