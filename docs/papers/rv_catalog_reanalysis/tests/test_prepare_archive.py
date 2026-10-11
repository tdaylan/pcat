import importlib.util
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "prepare_archive.py"
SPEC = importlib.util.spec_from_file_location("pcat_rv_archive_prep", SCRIPT_PATH)
PREPARATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARATION)


def test_radial_curve_parser_preserves_provenance_and_invalid_rows(tmp_path):
    curve_path = tmp_path / "UID_0000123_RVC_001.tbl"
    curve_path.write_text(
        '\\STAR_ID = "HD 123"\n'
        '\\NUMBER_OF_POINTS = "3"\n'
        '\\TIME_REFERENCE_FRAME = "BJD-TDB"\n'
        '\\BIBCODE = "2010ApJ...000..001A"\n'
        '\\INSTRUMENT = "HIRES"\n'
        "| BJD | RV | RV uncertainty |\n"
        "2450000.0 3.2 1.0\n"
        "2450001.0 NaN 1.2\n"
        "2450002.0 2.1 -1.0\n",
        encoding="utf-8",
    )

    curve, rows = PREPARATION._read_table_curve(curve_path, "RADIAL")

    assert curve["archive_uid"] == "0000123"
    assert curve["target_id"] == "HD 123"
    assert curve["time_reference_frame"] == "BJD-TDB"
    assert curve["bibcode"] == "2010ApJ...000..001A"
    assert len(rows) == 3
    assert curve["parse_errors"] == 0
    assert curve["point_count_matches_header"] is True
    assert PREPARATION._quality_reason(rows[0][:3]) == ""
    assert PREPARATION._quality_reason(rows[1][:3]) == "nonfinite_time_velocity_or_uncertainty"
    assert PREPARATION._quality_reason(rows[2][:3]) == "nonpositive_uncertainty"


def test_fulton_parser_preserves_per_row_instrument_codes(tmp_path):
    curve_path = tmp_path / "123456_rv.csv"
    curve_path.write_text(
        "# star HD number, 123456\n"
        "# k, pre-upgrade spectrograph\n"
        "hjd-2440000,rv,rv_err,inst\n"
        "10000.0,2.0,1.0,k\n"
        "10001.0,3.0,1.2,j\n",
        encoding="utf-8",
    )

    curve, rows = PREPARATION._read_fulton_curve(curve_path)

    assert curve["target_id"] == "HD 123456"
    assert curve["time_reference_frame"] == "HJD-2440000"
    assert [row[3] for row in rows] == ["k", "j"]
    assert curve["parse_errors"] == 0
    assert all(PREPARATION._quality_reason(row[:3]) == "" for row in rows)


def test_prepare_dataset_writes_curve_level_audit_and_observations(tmp_path):
    radial = tmp_path / "RADIAL"
    marvels = tmp_path / "MARVELS"
    fulton = tmp_path / "Fulton"
    output = tmp_path / "processed"
    radial.mkdir()
    marvels.mkdir()
    fulton.mkdir()
    (radial / "UID_0000123_RVC_001.tbl").write_text(
        '\\STAR_ID = "HD 123"\n\\NUMBER_OF_POINTS = "2"\n'
        '\\TIME_REFERENCE_FRAME = "JD"\n\\INSTRUMENT = "HIRES"\n'
        "2450000.0 1.0 1.0\n2450001.0 2.0 1.0\n",
        encoding="utf-8",
    )
    audit = PREPARATION.prepare_dataset(radial, marvels, fulton, output, {"test": True})

    assert audit["curve_counts"]["RADIAL"] == 1
    assert audit["observation_counts"]["RADIAL"]["valid_observations"] == 2
    curve_rows = list(__import__("csv").DictReader((output / "rv_curves.csv").open()))
    observation_rows = list(__import__("csv").DictReader((output / "rv_observations.csv").open()))
    assert curve_rows[0]["curve_id"] == "RADIAL:UID_0000123_RVC_001.tbl"
    assert len(observation_rows) == 2
    assert all(row["valid"] == "True" for row in observation_rows)
