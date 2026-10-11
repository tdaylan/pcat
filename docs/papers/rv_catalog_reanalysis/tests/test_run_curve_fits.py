import importlib.util
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "run_curve_fits.py"
SPEC = importlib.util.spec_from_file_location("pcat_rv_curve_fits", SCRIPT_PATH)
FITTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FITTER)


def _row(time, velocity, error, instrument="HIRES"):
    return {
        "time_value": str(time),
        "rv_mps": str(velocity),
        "rv_error_mps": str(error),
        "instrument": instrument,
    }


def test_identical_same_instrument_epoch_rows_are_removed_once():
    rows, removed, reason = FITTER._deduplicate_curve_rows([
        _row(2450000.0, 2.0, 1.0),
        _row(2450000.0, 2.0, 1.0),
        _row(2450001.0, 3.0, 1.1),
    ])

    assert len(rows) == 2
    assert removed == 1
    assert reason == ""


def test_conflicting_same_epoch_rows_are_not_combined():
    rows, removed, reason = FITTER._deduplicate_curve_rows([
        _row(2450000.0, 2.0, 1.0),
        _row(2450000.0, 2.3, 1.0),
    ])

    assert rows == []
    assert removed == 0
    assert reason == "conflicting_measurements_at_same_epoch"


def test_curve_run_names_and_analysis_directories_separate_seeds_and_settings(tmp_path):
    first_name = FITTER._run_name("RADIAL:curve.tbl", 17)
    second_name = FITTER._run_name("RADIAL:curve.tbl", 18)
    first_dir = FITTER._analysis_dir(tmp_path, 2000, 5, 10, 17)
    longer_dir = FITTER._analysis_dir(tmp_path, 100000, 5, 10, 17)
    second_seed_dir = FITTER._analysis_dir(tmp_path, 2000, 5, 10, 18)

    assert first_name != second_name
    assert first_dir != longer_dir
    assert first_dir != second_seed_dir


def test_empty_curve_is_skipped_before_importing_or_starting_pcat(tmp_path):
    curve = {
        "collection": "RADIAL",
        "curve_id": "RADIAL:empty",
        "target_id": "HD 123",
        "bibcode": "",
        "time_reference_frame": "JD",
    }

    result = FITTER._summarize_curve(curve, [], tmp_path, 2000, 5, 10)

    assert result["status"] == "skipped_no_valid_rows"
    assert result["reason"] == "no finite RV values with positive uncertainties"
    assert not (tmp_path / "analysis").exists()
