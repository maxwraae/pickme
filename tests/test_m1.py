import pytest
from pathlib import Path
from tn5verify.refs import load_folder
from tn5verify.grouping import build_groups
from tn5verify.characterize import annotate

FIXTURES = Path(__file__).parent / "fixtures" / "refs"

def test_load_refs():
    refs = load_folder(FIXTURES)
    names = {r.name for r in refs}
    assert names == {"rev340", "rev341", "rev342", "dcas9_v1", "empty_vec"}

def test_grouping_produces_correct_groups():
    refs = load_folder(FIXTURES)
    groups = build_groups(refs)
    # Expect: one 3-member group (rev family) + two singletons
    assert len(groups) == 3
    sizes = sorted([len(g.members) for g in groups])
    assert sizes == [1, 1, 3]
    rev_group = next(g for g in groups if len(g.members) == 3)
    assert set(rev_group.members) == {"rev340", "rev341", "rev342"}

def test_characterize_rev_family_has_variant_regions():
    refs = load_folder(FIXTURES)
    groups = build_groups(refs)
    groups = annotate(groups, refs)
    rev_group = next(g for g in groups if len(g.members) == 3)
    assert len(rev_group.variant_regions) >= 1

def test_characterize_singleton_has_no_variant_regions():
    refs = load_folder(FIXTURES)
    groups = build_groups(refs)
    groups = annotate(groups, refs)
    singletons = [g for g in groups if len(g.members) == 1]
    for s in singletons:
        assert len(s.variant_regions) == 0
        ref_name = s.members[0]
        # Singleton backbone spans whole sequence
        assert ref_name in s.backbone_intervals
        spans = s.backbone_intervals[ref_name]
        assert len(spans) >= 1

def test_features_parsed_from_gb():
    refs = load_folder(FIXTURES)
    rev340 = next(r for r in refs if r.name == "rev340")
    labels = {f.label for f in rev340.features}
    assert "AmpR" in labels
    assert "AmpR promoter" in labels
