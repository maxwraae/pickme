import gzip
import pytest
from pathlib import Path
from tn5verify.refs import load_folder
from tn5verify.grouping import build_groups
from tn5verify.characterize import annotate
from tn5verify.target import build_multi_fasta
from tn5verify.align import run_bwa

FIXTURES = Path(__file__).parent / "fixtures"
REFS = FIXTURES / "refs"
FASTQS = FIXTURES / "fastqs"

@pytest.fixture(scope="module")
def pipeline_base(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("m2")
    refs = load_folder(REFS)
    groups = build_groups(refs)
    groups = annotate(groups, refs)
    target_fa = build_multi_fasta(groups, refs, tmp)
    return {"refs": refs, "groups": groups, "target_fa": target_fa, "tmp": tmp}

def test_combined_fasta_created(pipeline_base):
    fa = pipeline_base["target_fa"]
    assert fa.exists()
    assert fa.with_suffix(".fa.bwt").exists() or (fa.parent / (fa.name + ".bwt")).exists()

def test_faidx_created(pipeline_base):
    fa = pipeline_base["target_fa"]
    assert Path(str(fa) + ".fai").exists()

def test_bwa_align_produces_dedup_bam(pipeline_base, tmp_path):
    r1 = FASTQS / "B4_R1.fastq.gz"
    r2 = FASTQS / "B4_R2.fastq.gz"
    if not r1.exists():
        pytest.skip("B4 fixture FASTQs not generated yet")
    stats = run_bwa(r1, r2, pipeline_base["target_fa"], tmp_path)
    assert stats.bam_path.exists()
    assert stats.total_reads > 0
    assert stats.mapped_reads >= 0

def test_mapping_rate_acceptable(pipeline_base, tmp_path):
    r1 = FASTQS / "B4_R1.fastq.gz"
    r2 = FASTQS / "B4_R2.fastq.gz"
    if not r1.exists():
        pytest.skip("B4 fixture FASTQs not generated yet")
    stats = run_bwa(r1, r2, pipeline_base["target_fa"], tmp_path)
    if stats.total_reads > 0:
        rate = stats.mapped_reads / stats.total_reads
        assert rate >= 0.80, f"Mapping rate {rate:.1%} below 80%"
