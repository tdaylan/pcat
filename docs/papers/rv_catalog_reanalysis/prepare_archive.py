"""Download and audit NASA Exoplanet Archive public RV curve collections.

The outputs preserve each archive curve as a separate series. They do not
combine curves from different references or silently resolve conflicting epochs.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import tarfile
import urllib.request

RADIAL_MANIFEST_URL = (
    "https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/wget_RADIAL.bat"
)
MARVELS_ARCHIVE_URL = (
    "https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/"
    "MARVELS_RV_bulk_download.tar.gz"
)
FULTON_ARCHIVE_URL = (
    "https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/"
    "Fulton_rv_data.tar.bz2"
)

OBSERVATION_FIELDS = (
    "collection", "curve_id", "archive_uid", "target_id", "target_alias",
    "bibcode", "reference", "telescope", "instrument", "observatory",
    "time_reference_frame", "time_value", "rv_mps", "rv_error_mps",
    "valid", "exclusion_reason",
)
CURVE_FIELDS = (
    "collection", "curve_id", "archive_uid", "target_id", "target_alias",
    "bibcode", "reference", "telescope", "instrument", "observatory",
    "time_reference_frame", "declared_observations", "parsed_observations",
    "valid_observations", "excluded_observations", "parse_errors",
    "duplicate_epochs_in_curve", "point_count_matches_header", "source_sha256",
)


def _download(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url, headers={"User-Agent": "PCAT-RV-archive-reanalysis/0.1"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        data = response.read()
    if not data:
        raise ValueError(f"Empty response from {url}")
    path.write_bytes(data)
    print(f"Downloaded {url} ({len(data)} bytes) to {path}")


def _extract_archive(archive_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved_root = output_dir.resolve()
    with tarfile.open(archive_path, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            output_path = (output_dir / member.name).resolve()
            if not output_path.is_relative_to(resolved_root):
                raise ValueError(f"Unsafe archive member path: {member.name}")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            input_file = archive.extractfile(member)
            if input_file is None:
                continue
            with input_file, output_path.open("wb") as stream:
                shutil.copyfileobj(input_file, stream)


def _ensure_sources(data_dir: Path) -> tuple[Path, Path, Path, str]:
    raw_dir = data_dir / "raw"
    radial_dir = raw_dir / "RADIAL"
    marvels_dir = raw_dir / "MARVELS"
    fulton_dir = raw_dir / "Fulton"
    manifest_path = raw_dir / "wget_RADIAL.bat"
    marvels_tar = raw_dir / "MARVELS_RV_bulk_download.tar.gz"
    fulton_tar = raw_dir / "Fulton_rv_data.tar.bz2"

    if not manifest_path.is_file():
        _download(RADIAL_MANIFEST_URL, manifest_path)
    manifest_text = manifest_path.read_text(encoding="utf-8", errors="replace")
    entries = {}
    for line in manifest_text.splitlines():
        match = re.match(r"wget\s+-O\s+(\S+)\s+(\S+)", line.strip())
        if match:
            filename, url = match.groups()
            entries[filename] = url.replace("http://", "https://").replace(":80/", "/")
    if not entries:
        raise ValueError(f"No RV curve downloads found in {manifest_path}")

    radial_dir.mkdir(parents=True, exist_ok=True)
    missing = [(filename, url) for filename, url in entries.items()
               if not (radial_dir / filename).is_file()]
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(_download, url, radial_dir / filename): filename
            for filename, url in missing
        }
        for future in as_completed(futures):
            future.result()

    if not marvels_tar.is_file():
        _download(MARVELS_ARCHIVE_URL, marvels_tar)
    if not fulton_tar.is_file():
        _download(FULTON_ARCHIVE_URL, fulton_tar)
    if not any(marvels_dir.rglob("*_RV.tbl")):
        _extract_archive(marvels_tar, marvels_dir)
    if not any(fulton_dir.rglob("*_rv.csv")):
        _extract_archive(fulton_tar, fulton_dir)

    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    return radial_dir, marvels_dir, fulton_dir, manifest_hash


def _read_table_curve(path: Path, collection: str) -> tuple[dict, list[tuple]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    metadata = {
        key: value.strip().strip('"').strip("'")
        for key, value in re.findall(r"^\\([A-Z_]+)\s*=\s*(.*)$", text, re.MULTILINE)
    }
    match = re.match(r"UID_(\d+)_RVC_(\d+)\.tbl$", path.name)
    archive_uid = match.group(1) if collection == "RADIAL" and match else ""
    curve_id = f"{collection}:{path.name}"
    rows = []
    parse_errors = 0
    for line in text.splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            time_value, velocity, uncertainty = map(float, fields[:3])
        except ValueError:
            if re.match(r"^\s*[+-]?(?:\d|\.\d)", line):
                parse_errors += 1
            continue
        rows.append((time_value, velocity, uncertainty, ""))
    if collection == "MARVELS":
        try:
            declared = int(metadata.get("NUMBER_OF_POINTS", ""))
        except ValueError:
            declared = None
        time_frame = metadata.get("TIME_REFERENCE_FRAME", "UNKNOWN")
        instrument = metadata.get("INSTRUMENT", "UNKNOWN")
        target_id = metadata.get("STAR_ID", path.stem)
        target_alias = path.stem
        bibcode = metadata.get("BIBCODE", "")
        reference = metadata.get("REFERENCE", "")
        telescope = metadata.get("TELESCOPE", "")
        observatory = metadata.get("OBSERVATORY_SITE", "")
    else:
        declared = None
        time_frame = metadata.get("TIME_REFERENCE_FRAME", "UNKNOWN")
        instrument = metadata.get("INSTRUMENT", "UNKNOWN")
        target_id = metadata.get("STAR_ID", archive_uid or path.stem)
        target_alias = metadata.get("STAR_ID", "")
        bibcode = metadata.get("BIBCODE", "")
        reference = metadata.get("REFERENCE", "")
        telescope = metadata.get("TELESCOPE", "")
        observatory = metadata.get("OBSERVATORY_SITE", "")
    valid_count = sum(
        all(map(__import__("math").isfinite, row[:3])) and row[2] > 0.0 for row in rows
    )
    duplicate_count = len(rows) - len({row[0] for row in rows})
    curve = {
        "collection": collection,
        "curve_id": curve_id,
        "archive_uid": archive_uid,
        "target_id": target_id,
        "target_alias": target_alias,
        "bibcode": bibcode,
        "reference": reference,
        "telescope": telescope,
        "instrument": instrument,
        "observatory": observatory,
        "time_reference_frame": time_frame,
        "declared_observations": declared,
        "parsed_observations": len(rows),
        "valid_observations": valid_count,
        "excluded_observations": len(rows) - valid_count,
        "parse_errors": parse_errors,
        "duplicate_epochs_in_curve": duplicate_count,
        "point_count_matches_header": declared is None or declared == len(rows),
        "source_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }
    return curve, rows


def _read_fulton_curve(path: Path) -> tuple[dict, list[tuple]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    target_match = re.search(r"^#\s*star HD number,\s*(\S+)\s*$", text, re.MULTILINE)
    target_id = f"HD {target_match.group(1)}" if target_match else path.stem
    descriptions = dict(re.findall(r"^#\s*([a-z])\s*,\s*(.+)$", text, re.MULTILINE))
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    rows = []
    parse_errors = 0
    for line in lines[1:]:
        fields = [field.strip() for field in line.split(",")]
        if len(fields) < 4:
            parse_errors += 1
            continue
        try:
            time_value, velocity, uncertainty = map(float, fields[:3])
        except ValueError:
            parse_errors += 1
            continue
        instrument_code = fields[3]
        rows.append((time_value, velocity, uncertainty, instrument_code))
    valid_count = sum(
        all(map(__import__("math").isfinite, row[:3])) and row[2] > 0.0 for row in rows
    )
    curve = {
        "collection": "FULTON",
        "curve_id": f"FULTON:{path.name}",
        "archive_uid": "",
        "target_id": target_id,
        "target_alias": path.stem,
        "bibcode": "2016PASP..128k4401H",
        "reference": "Howard and Fulton (2016)",
        "telescope": "Keck HIRES / Hamilton Spectrograph",
        "instrument": "; ".join(f"{key}={value}" for key, value in sorted(descriptions.items())),
        "observatory": "",
        "time_reference_frame": "HJD-2440000",
        "declared_observations": None,
        "parsed_observations": len(rows),
        "valid_observations": valid_count,
        "excluded_observations": len(rows) - valid_count,
        "parse_errors": parse_errors,
        "duplicate_epochs_in_curve": len(rows) - len({row[0] for row in rows}),
        "point_count_matches_header": True,
        "source_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }
    return curve, rows


def _quality_reason(row: tuple) -> str:
    import math

    if not all(math.isfinite(value) for value in row[:3]):
        return "nonfinite_time_velocity_or_uncertainty"
    if row[2] <= 0.0:
        return "nonpositive_uncertainty"
    return ""


def prepare_dataset(radial_dir: Path, marvels_dir: Path, fulton_dir: Path,
                    output_dir: Path, manifest_sha256: str) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "archive": "NASA Exoplanet Archive contributed RV time series",
        "archive_listing_last_updated": "2026-10-08",
        "radial_manifest_sha256": manifest_sha256,
        "curve_counts": {},
        "observation_counts": {},
    }
    epoch_curves = defaultdict(set)
    target_ids = defaultdict(set)
    time_frames = defaultdict(Counter)
    instrument_names = defaultdict(set)
    with (output_dir / "rv_curves.csv").open("w", newline="", encoding="utf-8") as curve_stream, \
         (output_dir / "rv_observations.csv").open("w", newline="", encoding="utf-8") as obs_stream:
        curve_writer = csv.DictWriter(curve_stream, fieldnames=CURVE_FIELDS)
        obs_writer = csv.DictWriter(obs_stream, fieldnames=OBSERVATION_FIELDS)
        curve_writer.writeheader()
        obs_writer.writeheader()
        source_groups = (
            ("RADIAL", sorted(radial_dir.rglob("UID_*_RVC_*.tbl"))),
            ("MARVELS", sorted(marvels_dir.rglob("*_RV.tbl"))),
            ("FULTON", sorted(fulton_dir.rglob("*_rv.csv"))),
        )
        for collection, paths in source_groups:
            counts = Counter()
            seen_curve_ids = set()
            for path in paths:
                if collection == "FULTON":
                    curve, rows = _read_fulton_curve(path)
                else:
                    curve, rows = _read_table_curve(path, collection)
                if curve["curve_id"] in seen_curve_ids:
                    continue
                seen_curve_ids.add(curve["curve_id"])
                curve_writer.writerow(curve)
                counts["curves"] += 1
                target_ids[collection].add(curve["target_id"])
                time_frames[collection][curve["time_reference_frame"]] += 1
                instrument_names[collection].add(curve["instrument"])
                counts["parsed_observations"] += len(rows)
                counts["valid_observations"] += curve["valid_observations"]
                counts["excluded_observations"] += curve["excluded_observations"]
                counts["parse_errors"] += curve["parse_errors"]
                counts["point_count_header_mismatches"] += int(
                    not curve["point_count_matches_header"]
                )
                for time_value, velocity, uncertainty, instrument_code in rows:
                    reason = _quality_reason((time_value, velocity, uncertainty))
                    instrument = instrument_code or curve["instrument"]
                    obs_writer.writerow({
                        "collection": collection,
                        "curve_id": curve["curve_id"],
                        "archive_uid": curve["archive_uid"],
                        "target_id": curve["target_id"],
                        "target_alias": curve["target_alias"],
                        "bibcode": curve["bibcode"],
                        "reference": curve["reference"],
                        "telescope": curve["telescope"],
                        "instrument": instrument,
                        "observatory": curve["observatory"],
                        "time_reference_frame": curve["time_reference_frame"],
                        "time_value": time_value,
                        "rv_mps": velocity,
                        "rv_error_mps": uncertainty,
                        "valid": not bool(reason),
                        "exclusion_reason": reason,
                    })
                    epoch_curves[(collection, curve["target_id"], curve["time_reference_frame"], round(time_value, 5))].add(curve["curve_id"])
            summary["curve_counts"][collection] = counts["curves"]
            summary["observation_counts"][collection] = dict(counts)
    summary["repeated_epoch_groups_across_curve_files"] = sum(
        len(curves) > 1 for curves in epoch_curves.values()
    )
    summary["repeated_epoch_extra_curve_records"] = sum(
        len(curves) - 1 for curves in epoch_curves.values() if len(curves) > 1
    )
    overlap_counts = Counter()
    for (collection, _, _, _), curves in epoch_curves.items():
        if len(curves) > 1:
            overlap_counts[collection] += 1
    summary["collection_metadata"] = {
        collection: {
            "unique_target_identifiers": len(target_ids[collection]),
            "curve_time_reference_counts": dict(time_frames[collection]),
            "instrument_labels": len(instrument_names[collection]),
            "repeated_epoch_groups_across_curves": overlap_counts[collection],
        }
        for collection in summary["curve_counts"]
    }
    summary["interpretation_note"] = (
        "Curves are kept separate. Overlapping measurements across publications are not combined or"
        " treated as independent system-level evidence by this preparation step."
    )
    (output_dir / "audit_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", type=Path, required=True,
        help="External storage directory for downloaded raw and normalized RV data.",
    )
    arguments = parser.parse_args()
    data_dir = arguments.data_dir.expanduser().resolve()
    repository_root = Path(__file__).resolve().parents[3]
    if data_dir.is_relative_to(repository_root):
        parser.error("--data-dir must be outside the repository to avoid storing raw RV archives in Git")
    radial_dir, marvels_dir, fulton_dir, manifest_sha256 = _ensure_sources(data_dir)
    summary = prepare_dataset(
        radial_dir, marvels_dir, fulton_dir, data_dir / "processed", manifest_sha256
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
