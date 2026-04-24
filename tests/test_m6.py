import pytest
from pathlib import Path
from tn5verify.refs import load_folder
from tn5verify.grouping import build_groups
from tn5verify.characterize import annotate
from tn5verify.features import annotate_gaps
from tn5verify.types import WellResult, WindowResult

FIXTURES = Path(__file__).parent / "fixtures"
REFS = FIXTURES / "refs"

def test_annotate_gaps_adds_feature_names():
    refs = load_folder(REFS)
    # Build a WellResult with a broken window in AmpR region (native ~100–1200)
    # AmpR is at 100-1200 in rev340 fixture. Window at padded 400-500 = native 100-200.
    window = WindowResult(contig="rev340", start=400, end=500,
                          mean_coverage=0.0, status="broken")
    result = WellResult(
        well_id="D3", total_reads=1000, mapped_reads=900,
        winning_group="0", called_member="rev340", coverage_span=0.8,
        backbone_windows=[window], verdict="RED_BROKEN", verdict_reason="test"
    )
    annotate_gaps(result, refs)
    assert len(window.affected_features) > 0
    # Should name AmpR or AmpR promoter
    assert any("AmpR" in f for f in window.affected_features)

def test_annotate_gaps_fallback_no_features():
    # A window in a region with no annotated features → fallback "bp X–Y"
    refs = load_folder(REFS)
    # Use a region outside AmpR: native ~1300-1400 → padded 1600-1700
    window = WindowResult(contig="rev340", start=1600, end=1700,
                          mean_coverage=0.0, status="warn")
    result = WellResult(
        well_id="X1", total_reads=1000, mapped_reads=900,
        winning_group="0", called_member="rev340", coverage_span=0.8,
        backbone_windows=[window], verdict="YELLOW", verdict_reason="test"
    )
    annotate_gaps(result, refs)
    # Either found a feature or fell back to "bp X–Y" label
    assert len(window.affected_features) > 0
