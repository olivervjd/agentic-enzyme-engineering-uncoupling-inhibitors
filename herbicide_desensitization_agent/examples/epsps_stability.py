"""Execute DDGun3D on the experimental Arabidopsis EPSPS apo structure.

This intentionally uses a documented custom ColabFold-MSA profile protocol,
not the upstream HHblits-search pipeline. Native DDGun scoring functions are
called without changing their equations, potentials or coefficients.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import importlib.util
import io
import json
import math
import re
import subprocess
import sys
from pathlib import Path

from herbicide_desensitization_agent.app.orchestrator.execution import ExecutionJournal, hash_artifacts

DDGUN_COMMIT = "bcf112c616bcd4eb9b18f193ec08cc704eb761b5"
KNOWN_CONTROLS = ("G177A", "T178I", "P182S")
HISTORICAL_HYPOTHESES = ("M288L", "M288I", "A329V", "A329S", "L356I", "L356M")
CANONICAL_OFFSET = 76
MATURE_LENGTH = 444
POTENTIALS = ("KYTJ820101", "HENS920102", "SKOJ970101", "BASU010101")
AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_a3m(text):
    """Remove A3M insertion columns; keep query-aligned gaps and all input rows."""
    rows, name, sequence = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if name is not None:
                rows.append((name, "".join(sequence)))
            name, sequence = line[1:].split()[0], []
        elif name is None:
            raise ValueError("MSA sequence appears before a FASTA header")
        else:
            sequence.append("".join(c for c in line if not c.islower() and c != "."))
    if name is not None:
        rows.append((name, "".join(sequence)))
    if len(rows) < 2 or not rows[0][1] or "-" in rows[0][1] or any(c not in AA for c in rows[0][1]):
        raise ValueError("A complete standard-amino-acid query and homologous alignment rows are required")
    query = rows[0][1]
    if any(len(sequence) != len(query) or any(c not in AA | {"X", "B", "Z", "U", "-"} for c in sequence) for _, sequence in rows):
        raise ValueError("A3M aligned rows must share query width and supported residue characters")
    if len({seq for _, seq in rows}) < 2:
        raise ValueError("A repeated query is not an evolutionary profile")
    return rows


def map_mutations(sequence, mutations=KNOWN_CONTROLS + HISTORICAL_HYPOTHESES, *, offset=CANONICAL_OFFSET):
    if len(sequence) != MATURE_LENGTH:
        raise ValueError("EPSPS protocol requires the exact 444-residue mature region")
    mapped = []
    for mutation in mutations:
        match = re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])", mutation)
        if not match:
            raise ValueError("Expected a single canonical point substitution")
        old, position, new = match.groups()
        local = int(position) - offset
        if not 1 <= local <= len(sequence) or sequence[local - 1] != old or old == new:
            raise ValueError(f"Mutation does not match the canonical-to-structure mapping: {mutation}")
        scope = "known_control" if mutation in KNOWN_CONTROLS else "historical_exploratory_hypothesis" if mutation in HISTORICAL_HYPOTHESES else "user_supplied_not_nominated"
        mapped.append({"mutation": mutation, "normalized_mutation": f"{old}{local}{new}", "canonical_position": int(position),
                       "normalized_position": local, "scope": scope, "new_recommendation": False})
    if len({row["mutation"] for row in mapped}) != len(mapped):
        raise ValueError("Duplicate mutations are not independent predictions")
    return mapped


def normalize_structure(cif, output, expected_sequence, *, source_chain="A"):
    from Bio.PDB import MMCIFParser, PDBIO
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    from Bio.SeqUtils import seq1

    source = MMCIFParser(QUIET=True).get_structure("7PXY", str(cif))
    if len(source) != 1 or source_chain not in source[0]:
        raise ValueError("Expected one experimental model containing the specified chain")
    residues = [r for r in source[0][source_chain] if r.id[0] == " " and "CA" in r]
    observed = "".join(seq1(r.resname, undef_code="X") for r in residues)
    if observed != expected_sequence or len(residues) != MATURE_LENGTH:
        raise ValueError("7PXY observed protein must exactly match the supplied mature WT MSA query; no alignment-based guessing")
    clean, model, chain = Structure("7PXY_normalized"), Model(0), Chain("A")
    clean.add(model)
    model.add(chain)
    mapping = []
    for i, residue in enumerate(residues, 1):
        new = Residue((" ", i, " "), residue.resname, residue.segid)
        atoms = {}
        for atom in residue.get_unpacked_list():
            if atom.element in ("H", "D"):
                continue
            rank = (atom.occupancy or 0, atom.altloc in (" ", "A"))
            if atom.name not in atoms or rank > atoms[atom.name][0]:
                atoms[atom.name] = rank, atom
        if not {"N", "CA", "C", "O"} <= set(atoms):
            raise ValueError(f"Incomplete experimental backbone at normalized position {i}")
        for _, atom in atoms.values():
            atom = atom.copy()
            atom.detach_parent()
            atom.set_altloc(" ")
            atom.disordered_flag = 0
            new.add(atom)
        chain.add(new)
        mapping.append({"chain": "A", "normalized_position": i, "canonical_position": i + CANONICAL_OFFSET,
                        "source_author_position": residue.id[1], "source_insertion_code": residue.id[2].strip(), "residue": observed[i-1]})
    writer = PDBIO()
    writer.set_structure(clean)
    writer.save(str(output))
    return {"chain": "A", "residue_count": len(residues), "canonical_range": [77, 520], "normalization_offset": CANONICAL_OFFSET,
            "sequence_sha256": hashlib.sha256(observed.encode()).hexdigest(), "mapping": mapping,
            "context": "experimental open apo conformation; waters and crystallization heterogens excluded"}


def compatibility_profile_module(ddgun_repo, directory):
    """Copy only ali2prof.py, removing an unused obsolete import and rU mode."""
    source = ddgun_repo / "tools" / "ali2prof.py"
    original = source.read_text()
    obsolete = "try:\n        from Bio.SearchIO._legacy import NCBIStandalone\nexcept:\n        from Bio.Blast import NCBIStandalone\n"
    if original.count(obsolete) != 1 or original.count('open(alifile, "rU")') != 1 or original.count("NCBIStandalone") != 2:
        raise ValueError("Upstream compatibility patch context changed; manual review required")
    patched = original.replace(obsolete, "# Removed obsolete unused NCBIStandalone import; scoring/profile math unchanged.\n").replace('open(alifile, "rU")', 'open(alifile, "r")')
    directory.mkdir(exist_ok=True)
    path = directory / "ali2prof.py"
    path.write_text(patched)
    diff = directory / "ali2prof.compatibility.patch"
    diff.write_text("".join(difflib.unified_diff(original.splitlines(True), patched.splitlines(True), fromfile="upstream/tools/ali2prof.py", tofile="compatibility/ali2prof.py")))
    return {"module": str(path), "patch": str(diff), "original_sha256": hashlib.sha256(original.encode()).hexdigest(),
            "patched_sha256": hashlib.sha256(patched.encode()).hexdigest(),
            "changes": ["remove unused obsolete Bio.NCBIStandalone import", "rU text mode to r; universal newlines remain Python default"],
            "algorithm_changes": False}


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_predictions(text, mapping):
    wanted = {r["normalized_mutation"]: r for r in mapping}
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    result = []
    for row in reader:
        mutation = row.get("VARIANT")
        if mutation not in wanted:
            raise ValueError("DDGun returned an unexpected mutation")
        raw = float(row["T_DDG[3D]"])
        if not math.isfinite(raw):
            raise ValueError("Non-finite DDGun estimate")
        effect = row["STABILITY[3D]"]
        expected = "Increase" if raw > 0 else "Decrease" if raw < 0 else "Neutral"
        if effect != expected:
            raise ValueError("DDGun physical sign does not match its native stability classification")
        result.append({**wanted[mutation], "ddgun3d_unfolding_ddg_kcal_mol": raw,
                       "folding_ddg_kcal_mol": -raw if raw else 0.0,
                       "ddgun_native_stability_effect": effect, "uncertainty_interval": None,
                       "calibration_status": "UNCALIBRATED_CUSTOM_MSA_PROTOCOL", "replicates": 1,
                       "evidence_state": "Computational stability estimate; no validated uncertainty interval",
                       "native_output_fields": row})
    if len(result) != len(wanted) or {r["normalized_mutation"] for r in result} != set(wanted):
        raise ValueError("DDGun did not return exactly one finite prediction per requested substitution")
    return result


def verify_outputs(directory):
    """Verify returned artifact bytes and sign mapping; export plot-ready rows.

    This is measured postflight validation, not another folding/stability run.
    Source tool execution timestamps remain in the original execution.json.
    """
    from herbicide_desensitization_agent.app.orchestrator.execution import validate_journal
    directory = Path(directory).resolve()
    result_path, raw_path = directory / "stability_results.json", directory / "ddgun3d_native.tsv"
    data = json.loads(result_path.read_text())
    original_journal = json.loads((directory / "execution.json").read_text())
    validate_journal(original_journal)
    verification = ExecutionJournal(directory / "verification_execution.json", run_id="epsps-ddgun-returned-artifact-verification")
    summary_path, csv_path = directory / "stability_summary.json", directory / "stability_summary.csv"
    def verify():
        if data["environment"]["ddgun_git_commit"] != DDGUN_COMMIT or data["protocol"] != "custom_ColabFoldMSA_DDGun3D" or data["calibrated"] is not False:
            raise ValueError("Returned stability protocol metadata differs from the executed custom protocol")
        observed = {}
        for artifact in data["output_artifacts"]:
            path = directory / Path(artifact["path"]).name
            actual = hash_artifacts([path])[0]
            if actual["sha256"] != artifact["sha256"] or actual["size_bytes"] != artifact["size_bytes"]:
                raise ValueError("Returned DDGun output bytes differ from execution provenance")
            observed[path.name] = actual
        query = "".join(line.strip() for line in (directory / "wt.fasta").read_text().splitlines() if not line.startswith(">"))
        mapping = map_mutations(query)
        if hashlib.sha256(query.encode()).hexdigest() != data["structure"]["sequence_sha256"]:
            raise ValueError("Returned structure/MSA sequence provenance mismatch")
        if data["structure"]["residue_count"] != MATURE_LENGTH or data["structure"]["canonical_range"] != [77, 520]:
            raise ValueError("Returned stability model has incorrect canonical sequence scope")
        parsed = parse_predictions(raw_path.read_text(), mapping)
        published = {row["mutation"]: row for row in data["mutations"]}
        if set(published) != {r["mutation"] for r in parsed}:
            raise ValueError("Returned mutation summary differs from the actual DDGun table")
        for row in parsed:
            for field in ("normalized_mutation", "scope", "ddgun3d_unfolding_ddg_kcal_mol", "folding_ddg_kcal_mol", "ddgun_native_stability_effect"):
                if row[field] != published[row["mutation"]][field]:
                    raise ValueError("Reported DDGun value/sign/scope differs from the native tool table")
        table_sha = observed["ddgun3d_native.tsv"]["sha256"]
        structures_sha = observed["7PXY_A_normalized_1_444.pdb"]["sha256"]
        rows = [{"mutation": row["mutation"], "canonical_position": row["canonical_position"],
                 "normalized_mutation": row["normalized_mutation"], "category": row["scope"],
                 "folding_ddg_estimate_kcal_mol": row["folding_ddg_kcal_mol"],
                 "unfolding_ddg_native_kcal_mol": row["ddgun3d_unfolding_ddg_kcal_mol"],
                 "uncertainty_low": None, "uncertainty_high": None, "uncertainty_status": "not calibrated",
                 "source_type": "computational_prediction", "tool": "DDGun3D", "tool_commit": DDGUN_COMMIT,
                 "protocol": data["protocol"], "source_structure": "7PXY chain A", "conformation": "open_apo",
                 "independent_structures": 1, "mutation_specific_structures_generated": 0,
                 "new_nomination": False, "native_output_sha256": table_sha, "structure_sha256": structures_sha,
                 "provenance": {"native_output": {"artifact": str(raw_path), "sha256": table_sha},
                                "structure": {"artifact": str(directory / "7PXY_A_normalized_1_444.pdb"), "sha256": structures_sha},
                                "profile": observed["wt.colabfold_profile.hssp"], "original_execution_journal": str(directory / "execution.json")}}
                for row in parsed]
        summary = {"schema_version": "1.0", "verification_status": "ARTIFACT_BYTES_VALUES_AND_SIGN_VERIFIED", "rows": rows,
                   "units": {"folding_ddg_estimate_kcal_mol": "kcal/mol", "unfolding_ddg_native_kcal_mol": "kcal/mol"},
                   "legend": {"folding_ddg_estimate_kcal_mol": "Positive means destabilizing; exact negative of native DDGun unfolding DDG",
                              "known_control": "Mapped prior control substitutions scored individually; no claim of target-specific experimental stability validation",
                              "historical_exploratory_hypothesis": "Archived hypotheses supplied for assessment; not newly generated or recommended mutations",
                              "uncertainty": "No calibrated intervals available; do not render missing intervals as zero-error estimates",
                              "structure": "All nine values use the same experimentally determined open apo WT structure; nine mutant folds were not generated"},
                   "limitations": data["limitations"], "counts": {"known_controls": 3, "historical_hypotheses": 6, "independent_structures": 1}}
        summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False))
        columns = [key for key in rows[0] if key != "provenance"]
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row[key] for key in columns})
        return summary
    try:
        return verification.execute("stability_artifact_verification", verify,
            inputs={"purpose": "verify returned actual tool outputs and export plot-ready provenance; no model rerun"},
            input_artifacts=[result_path, directory / "execution.json", raw_path], output_artifacts=[summary_path, csv_path],
            evidence_status="AVAILABLE", tools=[{"name": "deterministic DDGun output and artifact verifier", "source": str(Path(__file__).resolve())}])
    finally:
        verification.finalize()


def run(ddgun_repo, cif, msa, output, *, journal=None):
    ddgun_repo, cif, msa, output = map(lambda p: Path(p).resolve(), (ddgun_repo, cif, msa, output))
    output.mkdir(parents=True, exist_ok=True)
    if (output / "stability_results.json").exists():
        raise FileExistsError("Fresh stability runs must not overwrite prior scientific results")
    own_journal = journal is None
    journal = journal or ExecutionJournal(output / "execution.json", run_id="epsps-ddgun3d-custom-colabfoldmsa")
    try:
        def preflight():
            commit = subprocess.run(["git", "-C", str(ddgun_repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True).stdout.strip()
            if commit != DDGUN_COMMIT:
                raise ValueError("DDGun source does not match the reviewed pinned revision")
            subprocess.run(["git", "-C", str(ddgun_repo), "diff", "--exit-code", "--", "ddgun_3d.py", "tools", "data"], check=True, capture_output=True)
            import Bio, numpy
            return {"ddgun_git_commit": commit, "python": sys.version.split()[0], "biopython": Bio.__version__, "numpy": numpy.__version__}
        environment = journal.execute("stability_preflight", preflight, inputs={"ddgun_repo": str(ddgun_repo), "required_commit": DDGUN_COMMIT},
                                      input_artifacts=[cif, msa, ddgun_repo / "ddgun_3d.py"])
        alignment, fasta, profile = output / "colabfold_query_aligned.fasta", output / "wt.fasta", output / "wt.colabfold_profile.hssp"
        def prepare_alignment():
            rows = read_a3m(msa.read_text())
            mapping = map_mutations(rows[0][1])
            fasta.write_text(">EPSPS_WT_77_520\n" + rows[0][1] + "\n")
            with alignment.open("w") as handle:
                for i, (name, seq) in enumerate(rows):
                    label = "EPSPS_WT_77_520" if i == 0 else f"homolog_{i}_{name}"
                    handle.write(f">{label}\n{seq}\n")
            return {"query": rows[0][1], "aligned_depth": len(rows), "unique_aligned_sequences": len({s for _, s in rows}),
                    "mutation_mapping": mapping, "source": "existing ColabFold MSA; no new homolog search and no HHblits claim"}
        alignment_data = journal.execute("stability_msa_preparation", prepare_alignment, inputs={"source_msa": str(msa)}, input_artifacts=[msa], output_artifacts=[alignment, fasta])
        normalized_pdb = output / "7PXY_A_normalized_1_444.pdb"
        structure = journal.execute("stability_structure_normalization", lambda: normalize_structure(cif, normalized_pdb, alignment_data["query"]),
                                    inputs={"source": "7PXY", "chain": "A", "canonical_offset": CANONICAL_OFFSET}, input_artifacts=[cif, fasta], output_artifacts=[normalized_pdb])
        patch = journal.execute("stability_profile_compatibility", lambda: compatibility_profile_module(ddgun_repo, output / "compatibility"),
                                inputs={"allowed_changes": ["unused legacy import", "rU to r"]}, input_artifacts=[ddgun_repo / "tools" / "ali2prof.py"],
                                output_artifacts=lambda value: [value["module"], value["patch"]])
        def build_profile():
            module = load_module(patch["module"], "ddgun_ali2prof_compatibility")
            module.getfasta_profile(str(fasta), str(alignment), str(profile))
            if not profile.is_file() or profile.stat().st_size == 0:
                raise RuntimeError("DDGun ali2prof did not produce a profile")
            return {"profile": str(profile), "profile_generator": "upstream ali2prof.getfasta_profile", "msa_depth": alignment_data["aligned_depth"]}
        journal.execute("stability_profile_generation", build_profile, inputs={"protocol": "custom_ColabFoldMSA_DDGun3D"}, input_artifacts=[fasta, alignment, Path(patch["module"])], output_artifacts=[profile])
        dssp = ddgun_repo / "utils" / "dssp" / "dssp-2.0.4-linux-amd64"
        dssp_output = output / "7PXY_A.dssp"
        def calculate_dssp():
            command = [str(dssp), str(normalized_pdb), str(dssp_output)]
            process = subprocess.run(command, capture_output=True, text=True, timeout=300)
            (output / "dssp.stdout.log").write_text(process.stdout)
            (output / "dssp.stderr.log").write_text(process.stderr)
            process.check_returncode()
            return {"command": command, "returncode": process.returncode, "dssp_version": "2.0.4-linux-amd64"}
        journal.execute("stability_dssp", calculate_dssp, inputs={"tool": str(dssp)}, input_artifacts=[dssp, normalized_pdb],
                        output_artifacts=[dssp_output, output / "dssp.stdout.log", output / "dssp.stderr.log"])
        raw_table, feature_file = output / "ddgun3d_native.tsv", output / "ddgun3d_features.json"
        def score():
            sys.path.insert(0, str(ddgun_repo))
            ddgun = load_module(ddgun_repo / "ddgun_3d.py", "epsps_ddgun3d_pinned")
            mapped = alignment_data["mutation_mapping"]
            mutations = {m["normalized_mutation"]: [m["normalized_mutation"]] for m in mapped}
            try:
                features, profile_support, contacts = ddgun.get_muts_score(str(normalized_pdb), "A", str(dssp_output), str(profile),
                    mutations, list(POTENTIALS), [0.0, 5.0], win=2, outdir=str(output), add_seq=False)
            except SystemExit as exc:
                raise RuntimeError("DDGun scoring stopped without completing requested mutations") from exc
            text = "".join(ddgun.print_data(str(normalized_pdb), "A", features, profile_support, contacts, 2))
            rows = parse_predictions(text, mapped)
            raw_table.write_text(text)
            feature_file.write_text(json.dumps({"features": features, "profile_support": profile_support, "contact_environment": contacts},
                                               indent=2, allow_nan=False, default=lambda value: value.item()))
            return rows
        rows = journal.execute("stability", score, inputs={"mutations": alignment_data["mutation_mapping"], "potentials": POTENTIALS, "environment_angstrom": [0.0, 5.0], "sequence_window": 2},
                                input_artifacts=[normalized_pdb, dssp_output, profile, ddgun_repo / "ddgun_3d.py", *sorted((ddgun_repo / "data").glob("aaindex*.pkl"))],
                                output_artifacts=[raw_table, feature_file], evidence_status="AVAILABLE")
        result = {"schema_version": "1.0", "method": "DDGun3D", "protocol": "custom_ColabFoldMSA_DDGun3D",
                  "environment": environment, "structure": structure, "msa": {k: v for k, v in alignment_data.items() if k not in ("query", "mutation_mapping")},
                  "compatibility_patch": patch, "mutations": rows, "unit": "kcal/mol", "calibrated": False,
                  "sign_convention": {"native": "DDGun unfolding DDG; positive means stabilizing according to native tool effect output",
                                      "reported_folding": "folding_ddg = -ddgun_unfolding_ddg; positive folding DDG means destabilizing",
                                      "conversion": "Exact sign convention reversal, not conversion to binding affinity"},
                  "limitations": ["Custom ColabFold profile replaces upstream HHblits profile; this protocol has not been calibrated for EPSPS",
                                  "One experimentally determined open apo structure; no conformational ensemble or calibrated uncertainty interval",
                                  "No conclusion about ligand affinity, catalytic activity, native-function retention or herbicide resistance",
                                  "Historical substitutions remain archived hypotheses, not new nominations; known controls are scored singly",
                                  "Native tool rounds reported DDG to 0.1 kcal/mol; raw feature vectors are preserved",
                                  "Crystallization heterogens/waters excluded; no catalytic metal is assumed"],
                  "references": ["https://github.com/biofold/ddgun", "https://doi.org/10.1186/s12859-019-2923-1", "https://doi.org/10.1093/nar/gkac325"],
                  "input_artifacts": hash_artifacts([cif, msa]), "output_artifacts": hash_artifacts([normalized_pdb, alignment, fasta, profile, dssp_output, raw_table, feature_file])}
        result_path = output / "stability_results.json"
        journal.execute("stability_result_export", lambda: result_path.write_text(json.dumps(result, indent=2, allow_nan=False)),
                        inputs=result, output_artifacts=[result_path], evidence_status="AVAILABLE")
        return result
    finally:
        if own_journal:
            journal.finalize()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ddgun-repo", required=True)
    parser.add_argument("--structure-cif", required=True)
    parser.add_argument("--msa", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = run(args.ddgun_repo, args.structure_cif, args.msa, args.output)
    print(json.dumps({"status": "executed", "mutations_scored": len(result["mutations"]), "calibrated": result["calibrated"], "output": args.output}))


if __name__ == "__main__":
    main()
