"""Prompt contracts for configured structured reasoning adapters."""

SCIENTIFIC_PLANNER_PROMPT = """
Plan a fixed-registry Arabidopsis target-desensitization analysis. Verify the
AGI–herbicide pairing, select its registered mechanism-specific workflow, list
required biological context, and stop if required evidence is absent. Do not
propose mutations or treat model outputs as facts.
""".strip()

CANDIDATE_REVIEW_PROMPT = """
Review a computational candidate packet. Separate curated facts, model
predictions, assumptions, and unresolved uncertainty. Check native-function,
cofactor, complex, membrane, and conservation risks. Never approve a candidate
for experimental work; return a review recommendation for a human decision.
Numerical scientific outputs must come from validated tools or experiments;
never invent, estimate, or modify measurements or deterministic gate results.
""".strip()

DOMAIN_JUDGE_PROMPT = """
Judge biological correctness, evidence grounding, fixed-scope compliance,
constraint awareness, tool appropriateness, uncertainty calibration, negative
design logic, native-function preservation, actionability, governance,
reproducibility, and clarity. Require deterministic validator results and flag
unsupported claims. Numerical scientific outputs must come from validated
tools or experiments; do not supply missing measurements or override gates.
""".strip()

