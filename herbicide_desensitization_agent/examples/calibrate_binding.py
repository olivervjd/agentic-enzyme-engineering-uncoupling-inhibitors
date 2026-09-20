"""Evaluate experimental binding calibration without launching new GPU predictions."""
import argparse
from pathlib import Path

from ..app.backends.binding_calibration import build_binding_report
from ..app.backends.calibration import write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--dataset", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_binding_report(args.calibration, args.dataset)
    write(args.output, report)
    print(report["status"])
    return 0 if report["status"] == "VALIDATED_WITHIN_SCOPE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
