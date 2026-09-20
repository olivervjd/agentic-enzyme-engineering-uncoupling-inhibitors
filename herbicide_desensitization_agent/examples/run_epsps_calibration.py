"""Run the control-first EPSPS calibration protocol; no mutation design."""
import argparse
from pathlib import Path

from ..app.backends.calibration import EPSPSCalibrationRunner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--boltz-cache", type=Path, required=True)
    parser.add_argument("--boltz", default="boltz")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--reuse-calibration", type=Path)
    args = parser.parse_args()
    runner = EPSPSCalibrationRunner(args.previous, args.output, args.boltz_cache, args.boltz, args.reuse_calibration)
    if args.analyze_only:
        runner.analyze()
    else:
        runner.run()


if __name__ == "__main__":
    main()
