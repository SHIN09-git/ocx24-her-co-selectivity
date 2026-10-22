#!/usr/bin/env python3
"""Verify the release, audit saved outputs, or replay the frozen analysis offline."""
from pathlib import Path
import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
PREFIXES = ("00_admin/", "01_sources/", "02_data/", "03_analysis/", "04_results/")
FROZEN = ("preregistered_analysis_plan.json", "split_manifest.json", "protocol_lock.json")

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")

def verify_manifest():
    records = json.loads((ROOT / "RELEASE_MANIFEST.json").read_text(encoding="utf-8"))["files"]
    for rec in records:
        path = ROOT / rec["path"]
        if not path.is_file() or digest(path) != rec["sha256"]:
            raise RuntimeError("Release integrity mismatch: " + rec["path"])
    return records

def check_versions():
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError("Full replay requires Python 3.12.")
    versions = {}
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        package, expected = line.split("==")
        actual = importlib.metadata.version(package)
        if actual != expected:
            raise RuntimeError(f"{package}: expected {expected}, found {actual}")
        versions[package] = actual
    return versions

def run_step(stage, script):
    print("Running " + script, flush=True)
    started = time.monotonic()
    log = stage / "logs" / (Path(script).stem + ".log")
    env = os.environ.copy()
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    with log.open("w", encoding="utf-8") as out:
        result = subprocess.run(
            [sys.executable, "-E", "-s", "-B", "-X", "utf8", str(stage / "03_analysis" / script)],
            cwd=stage / "03_analysis", env=env, stdout=out, stderr=subprocess.STDOUT,
        )
    seconds = round(time.monotonic() - started, 3)
    if result.returncode:
        tail = log.read_text(encoding="utf-8")[-6000:]
        raise RuntimeError(f"{script} failed (exit {result.returncode}). See {log}\n{tail}")
    print(f"Completed {script} ({seconds:.1f} s)", flush=True)
    return {"script": script, "seconds": seconds, "status": "passed"}

def review_summary(stage):
    result = json.loads((stage / "04_results/authorized_review/summary.json").read_text(encoding="utf-8"))
    expected = {
        "raw_pair_records_checked": 441, "designated_key_records_checked": 20,
        "mapping_anomalies_checked": 6, "batch_date_conflict_electrodes_checked": 216,
        "exclusions_checked": 1, "saved_metric_tasks_checked": 23, "models_refitted": 0,
        "frozen_inputs_and_predictions_unchanged": True,
        "status": "source_and_calculation_checks_passed",
    }
    for field, value in expected.items():
        if result[field] != value:
            raise RuntimeError(f"Unexpected audit count: {field} = {result[field]!r}")
    return result

class Comparison:
    def __init__(self):
        self.numeric_values = 0
        self.max_abs_difference = 0.0

    def equal(self, actual, expected, where):
        if isinstance(expected, dict):
            if not isinstance(actual, dict) or actual.keys() != expected.keys():
                raise RuntimeError("JSON keys differ: " + where)
            for key in expected:
                self.equal(actual[key], expected[key], where + "/" + key)
        elif isinstance(expected, list):
            if not isinstance(actual, list) or len(actual) != len(expected):
                raise RuntimeError("List length differs: " + where)
            for index, (a, b) in enumerate(zip(actual, expected)):
                self.equal(a, b, f"{where}/{index}")
        elif isinstance(expected, bool) or expected is None:
            if actual != expected:
                raise RuntimeError("Value differs: " + where)
        else:
            try:
                a, b = float(actual), float(expected)
            except (TypeError, ValueError):
                if actual != expected:
                    raise RuntimeError(f"Value differs: {where}: {actual!r} != {expected!r}")
                return
            self.numeric_values += 1
            if math.isnan(a) and math.isnan(b):
                return
            if not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-9):
                raise RuntimeError(f"Number differs: {where}: {a!r} != {b!r}")
            if math.isfinite(a) and math.isfinite(b):
                self.max_abs_difference = max(self.max_abs_difference, abs(a-b))

    def file(self, actual, expected):
        if expected.suffix == ".json":
            self.equal(json.loads(actual.read_text(encoding="utf-8")),
                       json.loads(expected.read_text(encoding="utf-8")), expected.name)
            return
        def read(path):
            with path.open(encoding="utf-8-sig", newline="") as stream:
                reader = csv.DictReader(stream)
                columns = reader.fieldnames
                rows = list(reader)
            # Aggregate task order can differ across filesystems.
            if "task" in (columns or []):
                rows.sort(key=lambda row: tuple(row.get(k, "") for k in
                          ("task", "fold", "model", "variant", "electrode_id")))
            return columns, rows
        a_columns, a_rows = read(actual)
        b_columns, b_rows = read(expected)
        if a_columns != b_columns:
            raise RuntimeError("CSV columns differ: " + expected.name)
        self.equal(a_rows, b_rows, expected.name)

def compare_replay(stage):
    comparison = Comparison()
    checked = []
    for folder in ("02_data", "04_results"):
        for original in sorted((ROOT / folder).rglob("*")):
            if original.suffix not in (".json", ".csv"):
                continue
            relative = original.relative_to(ROOT)
            actual = stage / relative
            if not actual.is_file():
                raise RuntimeError("Missing reproduced output: " + relative.as_posix())
            comparison.file(actual, original)
            checked.append(relative.as_posix())
    return {"files_compared": len(checked), "files": checked,
            "numeric_values_compared": comparison.numeric_values,
            "max_abs_difference": comparison.max_abs_difference,
            "absolute_tolerance": 1e-9, "relative_tolerance": 1e-12,
            "excluded": "Archived figure pixels; timestamped replay records; generated model binaries."}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify-only", action="store_true", help="Check release hashes; write nothing.")
    mode.add_argument("--full", action="store_true", help="Rebuild data and refit all frozen models.")
    parser.add_argument("--output", type=Path, help="New output directory (must not exist).")
    args = parser.parse_args()
    records = verify_manifest()
    print(f"Release integrity: {len(records)} files passed.", flush=True)
    if args.verify_only:
        return
    versions = check_versions() if args.full else {}
    core = {x["path"]: x["sha256"] for x in records if x["path"].startswith(PREFIXES)
            or x["path"] in ("reproduce.py", "requirements.txt")}
    core_sha = hashlib.sha256(json.dumps(core, sort_keys=True).encode()).hexdigest()
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    stage = (args.output or ROOT / "reproductions" / stamp).resolve()
    if stage == ROOT or (stage.is_relative_to(ROOT) and not stage.is_relative_to(ROOT / "reproductions")):
        raise RuntimeError("Output must be outside the release, or under reproductions/.")
    if stage.exists():
        raise RuntimeError("Output exists; select a new directory: " + str(stage))
    stage.mkdir(parents=True)
    (stage / "logs").mkdir()
    for folder in ("00_admin", "01_sources", "02_data", "03_analysis", "04_results"):
        if args.full and folder in ("00_admin", "02_data", "04_results"):
            (stage / folder).mkdir()
        else:
            shutil.copytree(ROOT / folder, stage / folder)
    report = {
        "mode": "full_replay" if args.full else "saved_output_audit",
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "python": platform.python_version(), "platform": platform.platform(),
        "dependencies": versions, "release_files_verified_at_start": len(records),
        "analysis_input_fingerprint_sha256": core_sha, "analysis_input_sha256": core,
        "scope": "Reproduction of a completed public-data analysis; not a new blind evaluation or peer review.",
        "steps": [], "status": "running",
    }
    started = time.monotonic()
    try:
        if args.full:
            report["steps"].append(run_step(stage, "build_dataset.py"))
            for name in ("features.csv", "development_targets.csv", "holdout_targets.csv"):
                if digest(stage / "02_data" / name) != digest(ROOT / "02_data" / name):
                    raise RuntimeError("Frozen data rebuild differs: " + name)
            split = json.loads((stage / "00_admin/split_manifest.json").read_text())
            expected = json.loads((ROOT / "00_admin/split_manifest.json").read_text())
            split.pop("created_utc"); expected.pop("created_utc")
            if split != expected:
                raise RuntimeError("Frozen split differs.")
            shutil.copy2(stage / "00_admin/split_manifest.json", stage / "rebuilt_split_manifest.json")
            for name in FROZEN:
                shutil.copy2(ROOT / "00_admin" / name, stage / "00_admin" / name)
            report["frozen_input_rebuild"] = "Three CSV files byte-identical; split and original lock preserved."
            for script in ("run_primary.py", "run_controls.py", "evaluate_holdout.py",
                           "validate_analysis.py", "supplemental_audit.py",
                           "audit_coverage_and_screening.py", "screening_uncertainty.py",
                           "export_family_table.py"):
                report["steps"].append(run_step(stage, script))
            report["comparison"] = compare_replay(stage)
        report["steps"].append(run_step(stage, "authorized_review.py"))
        report["source_and_metric_checks"] = review_summary(stage)
        for relative, expected in core.items():
            if digest(ROOT / relative) != expected:
                raise RuntimeError("Release input changed during replay: " + relative)
        report["release_inputs_unchanged"] = True
        report["status"] = "passed"
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
        raise
    finally:
        report["elapsed_seconds"] = round(time.monotonic() - started, 3)
        report["completed_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        save(stage / "reproduction_report.json", report)
    print(json.dumps({"status": report["status"], "mode": report["mode"],
                      "output": str(stage), "elapsed_seconds": report["elapsed_seconds"]},
                     ensure_ascii=False), flush=True)

if __name__ == "__main__":
    main()

