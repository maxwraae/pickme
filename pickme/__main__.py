"""CLI entry: `pickme` (TUI) and `pickme run <name>` (headless)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pickme.runs import Run, scan_input_dir


def _parse_insert_region(text: str) -> tuple[int, int]:
    try:
        s, e = text.split("-", 1)
        return int(s), int(e)
    except ValueError as exc:
        raise SystemExit(
            f"--insert-region must be 'START-END', got {text!r}"
        ) from exc


def _cmd_run(args: argparse.Namespace) -> int:
    input_root = Path(args.input_dir).resolve()
    output_root = Path(args.output_dir).resolve()

    matches = [
        r for r in scan_input_dir(input_root, output_root) if r.name == args.run_name
    ]
    if not matches:
        print(
            f"ERROR: no run named {args.run_name!r} under {input_root}",
            file=sys.stderr,
        )
        return 2

    run: Run = matches[0]
    if run.status == "invalid":
        print(
            f"ERROR: run {run.name!r} is invalid: {run.summary}",
            file=sys.stderr,
        )
        return 2

    insert_region = None
    if args.insert_region:
        insert_region = _parse_insert_region(args.insert_region)

    from pickme.runner import run_pipeline

    def _progress(well_id: str, stage: str, payload: dict) -> None:
        if stage == "start":
            print(f"Aligning {payload['n_wells']} wells...")
        elif stage == "align":
            print(
                f"  {well_id}: {payload.get('mapped_reads', 0)}"
                f"/{payload.get('total_reads', 0)} mapped"
            )
        elif stage == "identify":
            extra = ""
            if payload.get("flagged"):
                extra = f"  scan:{payload['flagged']} flagged"
            print(
                f"  {well_id}: {payload.get('verdict', '?')} "
                f"({payload.get('called_member', 'unknown')}){extra}"
            )
        elif stage == "summary":
            print("\nVerdict summary:")
            for v, c in sorted(payload.get("verdicts", {}).items()):
                print(f"  {v}: {c}")
            print(f"Outputs: {payload.get('xlsx')}")

    run_pipeline(
        run,
        scan_mode=args.scan_mode,
        insert_region=insert_region,
        threads=args.threads,
        on_progress=_progress,
    )
    return 0


def _cmd_tui(args: argparse.Namespace) -> int:
    from pickme.app import PickmeApp

    PickmeApp(
        input_dir=Path(args.input_dir),
        output_dir=Path(args.output_dir),
        scan_mode=args.scan_mode,
    ).run()
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="pickme",
        description="Tn5 construct verification for 96-well plates.",
    )
    parser.add_argument(
        "--input-dir",
        default="input",
        help="Root folder containing per-run subdirs (default: ./input).",
    )
    parser.add_argument(
        "--output-dir",
        default="output",
        help="Root folder for results (default: ./output).",
    )

    subparsers = parser.add_subparsers(dest="cmd")

    tui_p = subparsers.add_parser("tui", help="Launch the TUI (default).")
    tui_p.add_argument(
        "--scan-mode",
        choices=["quick", "insert", "whole"],
        default="quick",
    )

    run_p = subparsers.add_parser("run", help="Run headless on a single run name.")
    run_p.add_argument("run_name", help="Name of a subfolder under input/")
    run_p.add_argument(
        "--scan-mode",
        choices=["quick", "insert", "whole"],
        default="quick",
    )
    run_p.add_argument(
        "--insert-region",
        default=None,
        help="START-END, required when --scan-mode=insert.",
    )
    run_p.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Alignment threads per well (default: half of CPU count).",
    )

    args = parser.parse_args()

    if args.cmd == "run":
        sys.exit(_cmd_run(args))

    # Default: TUI (either `pickme` or `pickme tui`)
    if args.cmd is None:
        # argparse for tui subparser hasn't run; inject defaults
        args.scan_mode = "quick"
    sys.exit(_cmd_tui(args))


if __name__ == "__main__":
    main()
