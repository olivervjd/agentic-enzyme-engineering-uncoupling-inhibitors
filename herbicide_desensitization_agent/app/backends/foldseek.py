from __future__ import annotations

import shutil
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
                "query,target,qtmscore,ttmscore,lddt,alnlen,qlen",
            ], check=True, capture_output=True, text=True)
            line = result.read_text(encoding="utf-8").splitlines()[0].split("\t")
        return {
            "query_tm_score": float(line[2]), "target_tm_score": float(line[3]),
            "alignment_lddt": float(line[4]), "alignment_coverage": int(line[5]) / int(line[6]),
        }
