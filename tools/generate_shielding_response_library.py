#!/usr/bin/env python3
"""Generate a reusable WEP thickness/dose response table.

The normal dynamic MOUSE calculation writes ``shielding_inputs.json`` and
``shielding_leakage_source.json`` in its shielding run directory.  This tool
reuses those files to evaluate additional thicknesses without repeating the
detailed core or vessel-transfer calculations.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from core_design.shielding import evaluate_shielding_candidate


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate shutdown dose at specified shield thicknesses using a "
            "source exported by a previous dynamic MOUSE shielding run."
        )
    )
    parser.add_argument(
        "--inputs",
        required=True,
        type=Path,
        help="Path to shielding_inputs.json",
    )
    parser.add_argument(
        "--source",
        required=True,
        type=Path,
        help="Path to shielding_leakage_source.json",
    )
    parser.add_argument(
        "--thicknesses",
        required=True,
        type=float,
        nargs="+",
        help="Shield thicknesses to evaluate, in cm",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("shielding_response_library"),
        help="Directory for OpenMC runs and CSV/JSON output",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    params = json.loads(args.inputs.read_text(encoding="utf-8"))
    source_record = json.loads(args.source.read_text(encoding="utf-8"))
    leakage_source = source_record.get("outside_intake_source", source_record)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for thickness in sorted(set(args.thicknesses)):
        print(f"Evaluating {thickness:.3f} cm...")
        result = evaluate_shielding_candidate(
            params,
            leakage_source,
            thickness,
            output_dir / f"candidate_{thickness:.2f}_cm",
        )
        rows.append(result)
        print(
            f"  {result['dose_mrem_per_h']:.6g} +/- "
            f"{result['dose_std_dev_mrem_per_h']:.3g} mrem/h"
        )

    json_path = output_dir / "shielding_response_library.json"
    json_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    csv_path = output_dir / "shielding_response_library.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

