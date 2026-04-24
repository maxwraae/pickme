import pytest
from pathlib import Path
import tempfile
from tn5verify.refs import load_folder
from tn5verify.grouping import build_groups
from tn5verify.characterize import annotate
from tn5verify.target import build_multi_fasta
from tn5verify.align import run_bwa
from tn5verify.identify import classify

FIXTURES = Path(__file__).parent / "fixtures"
REFS = FIXTURES / "refs"
FASTQS = FIXTURES / "fastqs"

@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("m3")
    refs = load_folder(REFS)
    groups = build_groups(refs)
    groups = annotate(groups, refs)
    target_fa = build_multi_fasta(groups, refs, tmp)
    return {"refs": refs, "groups": groups, "target_fa": target_fa, "tmp": tmp}

def test_classify_rev340_well(pipeline, tmp_path):
    r1 = FASTQS / "B4_R1.fastq.gz"
    r2 = FASTQS / "B4_R2.fastq.gz"
    stats = run_bwa(r1, r2, pipeline["target_fa"], tmp_path)
    result = classify(stats.bam_path, pipeline["groups"], pipeline["refs"])
    assert result.called_member is not None
    assert "rev" in result.called_member  # should call one of the rev family

def test_classify_unknown_well(pipeline, tmp_path):
    r1 = FASTQS / "G5_R1.fastq.gz"
    r2 = FASTQS / "G5_R2.fastq.gz"
    if not r1.exists():
        pytest.skip("G5 fixture not generated")
    stats = run_bwa(r1, r2, pipeline["target_fa"], tmp_path)
    result = classify(stats.bam_path, pipeline["groups"], pipeline["refs"])
    # Random reads → NO_DATA or RED_UNKNOWN
    assert result.verdict in ("NO_DATA", "RED_UNKNOWN")

def test_classify_returns_well_result(pipeline, tmp_path):
    r1 = FASTQS / "B4_R1.fastq.gz"
    r2 = FASTQS / "B4_R2.fastq.gz"
    stats = run_bwa(r1, r2, pipeline["target_fa"], tmp_path)
    result = classify(stats.bam_path, pipeline["groups"], pipeline["refs"])
    assert result.total_reads >= 0
    assert result.mapped_reads >= 0
    assert result.coverage_span >= 0.0
    assert result.winning_group is not None
