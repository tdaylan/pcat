"""Run resumable, independent PCAT fits for prepared archive RV curves.

Curves are not combined by target in this first-pass runner. This avoids
silently multiplying overlapping published reductions. Results are descriptive
per-curve posteriors, not population-level occurrence estimates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np

RESULT_FIELDS = (
    "collection", "curve_id", "seed", "target_id", "bibcode", "time_reference_frame",
    "valid_observations", "duplicate_rows_removed", "status", "reason",
    "runtime_seconds", "draws",
    "map_planet_count", "p_n0", "p_n1", "p_n2", "p_n3", "p_n4", "p_n5",
)


def _run_name(curve_id: str, seed: int) -> str:
    digest = hashlib.sha256(f"{seed}:{curve_id}".encode("utf-8")).hexdigest()[:20]
    return f"rvcurve_{digest}"


def _analysis_dir(data_dir: Path, number_sweeps: int, maximum_planets: int,
                  minimum_observations: int, seed: int) -> Path:
    run_config = (
        f"seed_{seed}_sweeps_{number_sweeps}_maxn_{maximum_planets}"
        f"_minobs_{minimum_observations}"
    )
    return data_dir / "analysis" / run_config


def _deduplicate_curve_rows(rows: list[dict]) -> tuple[list[dict], int, str]:
    ordered_rows = sorted(rows, key=lambda row: float(row["time_value"]))
    unique_rows = []
    number_removed = 0
    index = 0
    while index < len(ordered_rows):
        first_row = ordered_rows[index]
        epoch = float(first_row["time_value"])
        same_epoch = [first_row]
        index += 1
        while (index < len(ordered_rows)
               and abs(float(ordered_rows[index]["time_value"]) - epoch) <= 1e-6):
            same_epoch.append(ordered_rows[index])
            index += 1
        fingerprints = {
            (row["instrument"], row["rv_mps"], row["rv_error_mps"])
            for row in same_epoch
        }
        if len(fingerprints) != 1:
            return [], number_removed, "conflicting_measurements_at_same_epoch"
        unique_rows.append(first_row)
        number_removed += len(same_epoch) - 1
    return unique_rows, number_removed, ""


def _summarize_curve(curve: dict, rows: list[dict], data_dir: Path,
                     number_sweeps: int, maximum_planets: int,
                     minimum_observations: int, seed: int = 17) -> dict:
    result = {key: "" for key in RESULT_FIELDS}
    for key in ("collection", "curve_id", "target_id", "bibcode", "time_reference_frame"):
        result[key] = curve.get(key, "")
    result["seed"] = seed
    valid_rows = [row for row in rows if row["valid"]]
    result["valid_observations"] = len(valid_rows)
    result["runtime_seconds"] = ""
    result["draws"] = 0
    result["duplicate_rows_removed"] = 0
    if not valid_rows:
        result["status"] = "skipped_no_valid_rows"
        result["reason"] = "no finite RV values with positive uncertainties"
        return result
    if len(valid_rows) < minimum_observations:
        result["status"] = "skipped_too_few_rows"
        result["reason"] = f"{len(valid_rows)} valid rows; minimum is {minimum_observations}"
        return result

    valid_rows, duplicates_removed, duplicate_error = _deduplicate_curve_rows(valid_rows)
    result["duplicate_rows_removed"] = duplicates_removed
    if duplicate_error:
        result["status"] = "skipped_conflicting_duplicate_epoch"
        result["reason"] = duplicate_error
        return result
    if len(valid_rows) < minimum_observations:
        result["status"] = "skipped_too_few_after_deduplication"
        result["reason"] = f"{len(valid_rows)} rows remain; minimum is {minimum_observations}"
        return result

    analysis_dir = _analysis_dir(
        data_dir, number_sweeps, maximum_planets, minimum_observations, seed,
    )
    run_name = _run_name(curve["curve_id"], seed)
    from pcat.main import readfile, retr_pathrun
    run_root = Path(retr_pathrun(str(analysis_dir), run_name))
    posterior_base = run_root / "data" / "outp" / run_name / "gdatfinlpost"
    if not posterior_base.with_suffix(".p").is_file() or not posterior_base.with_suffix(".h5").is_file():
        output_path = run_root / "data" / "outp" / run_name
        if output_path.exists():
            result["status"] = "incomplete_existing_run"
            result["reason"] = "run cache exists without a complete final posterior; preserved"
            return result

    time_values = np.asarray([float(row["time_value"]) for row in valid_rows])
    velocities = np.asarray([float(row["rv_mps"]) for row in valid_rows])
    uncertainties = np.asarray([float(row["rv_error_mps"]) for row in valid_rows])
    instruments = np.asarray([row["instrument"] or "UNKNOWN" for row in valid_rows])
    order = np.argsort(time_values)
    time_values = time_values[order]
    velocities = velocities[order]
    uncertainties = uncertainties[order]
    instruments = instruments[order]
    if np.any(np.diff(time_values) <= 1e-6):
        result["status"] = "skipped_duplicate_epoch"
        result["reason"] = "PCAT RV input requires distinct epochs; overlap preserved for separate curation"
        return result

    if posterior_base.with_suffix(".p").is_file() and posterior_base.with_suffix(".h5").is_file():
        posterior = readfile(str(posterior_base))
    else:
        from pcat import sampling
        from pcat.radial_velocity import retr_dictpcatrvel

        _, instrument_index = np.unique(instruments, return_inverse=True)
        burn_count = number_sweeps // 3
        configuration = retr_dictpcatrvel(
            time_values,
            velocities,
            uncertainties,
            instrument_index,
                str(analysis_dir),
            run_name,
            maxmnumbplan=maximum_planets,
            boollinetren=False,
            probtran=0.7,
            probspmr=0.0,
            probjump=0.2,
                typeseed=seed,
            inittype="rand",
            boolburntmpr=True,
            factburntmpr=0.8,
            stdvpropelemfire=[1e-2, 1e-4, 3e-2, 3e-2, 3e-2],
            numbswep=number_sweeps,
            numbburn=burn_count,
            numbsamp=max(number_sweeps // 50, 100),
            numbproc=1,
            boolmakeplot=False,
            boolmakeplotinit=False,
            makeanim=False,
            typeverb=0,
        )
        start = time.perf_counter()
        sampling.sample(**configuration)
        result["runtime_seconds"] = f"{time.perf_counter() - start:.6f}"
        if not posterior_base.with_suffix(".p").is_file() or not posterior_base.with_suffix(".h5").is_file():
            result["status"] = "failed_no_final_state"
            result["reason"] = "PCAT returned without a paired final posterior state"
            return result
        posterior = readfile(str(posterior_base))

    planet_counts = np.asarray(posterior.listpostnumbelem).astype(int).reshape(-1)
    frequencies = np.bincount(planet_counts, minlength=maximum_planets + 1)[:maximum_planets + 1]
    result["draws"] = int(planet_counts.size)
    result["map_planet_count"] = int(np.argmax(frequencies))
    for count, frequency in enumerate(frequencies):
        result[f"p_n{count}"] = f"{frequency / max(planet_counts.size, 1):.8f}"
    result["status"] = "completed"
    return result


def run_archive_fits(data_dir: Path, number_sweeps: int, maximum_planets: int,
                     minimum_observations: int, limit: int | None = None,
                     seed: int = 17) -> Path:
    processed_dir = data_dir / "processed"
    curve_path = processed_dir / "rv_curves.csv"
    observation_path = processed_dir / "rv_observations.csv"
    if not curve_path.is_file() or not observation_path.is_file():
        raise FileNotFoundError("Run prepare_archive.py first to create normalized observations")
    with curve_path.open(newline="", encoding="utf-8") as stream:
        curves = list(csv.DictReader(stream))
    grouped = {}
    with observation_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            row["valid"] = row["valid"].lower() == "true"
            grouped.setdefault(row["curve_id"], []).append(row)
    analysis_dir = _analysis_dir(
        data_dir, number_sweeps, maximum_planets, minimum_observations, seed,
    )
    result_path = analysis_dir / "per_curve_planet_count.csv"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    completed_ids = set()
    if result_path.is_file():
        with result_path.open(newline="", encoding="utf-8") as stream:
            completed_ids = {row["curve_id"] for row in csv.DictReader(stream)}
    pending = [curve for curve in curves if curve["curve_id"] not in completed_ids]
    if limit is not None:
        pending = pending[:limit]
    write_header = not result_path.is_file()
    with result_path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_FIELDS)
        if write_header:
            writer.writeheader()
        for index, curve in enumerate(pending, 1):
            result = _summarize_curve(
                curve, grouped.get(curve["curve_id"], []), data_dir,
                number_sweeps, maximum_planets, minimum_observations, seed,
            )
            writer.writerow(result)
            stream.flush()
            print(f"{index}/{len(pending)} {curve['curve_id']}: {result['status']} {result['reason']}")
    return result_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--numbswep", type=int, default=100_000)
    parser.add_argument("--maxmnumbplan", type=int, default=5)
    parser.add_argument("--min-observations", type=int, default=10)
    parser.add_argument("--limit", type=int, default=None,
                        help="Process the first N pending curves for a smoke/pilot run")
    parser.add_argument("--seed", type=int, default=17)
    arguments = parser.parse_args()
    output = run_archive_fits(
        arguments.data_dir.expanduser().resolve(), arguments.numbswep,
        arguments.maxmnumbplan, arguments.min_observations, arguments.limit,
        arguments.seed,
    )
    print(f"Writing to {output}...")


if __name__ == "__main__":
    main()
