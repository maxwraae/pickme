from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


@dataclass(frozen=True)
class Feature:
    kind: str        # "CDS", "promoter", "terminator", ...
    label: str       # "AmpR", "T7 promoter", ...
    start: int       # native coords, 0-based
    end: int         # exclusive


@dataclass(frozen=True)
class Reference:
    name: str
    sequence: str    # native, unpadded
    features: tuple[Feature, ...]   # use tuple not list (frozen dataclass)
    source_path: Path


@dataclass
class VariantRegion:
    # coordinates in EACH member's own native space
    member_spans: dict[str, tuple[int, int]]  # {ref_name: (start, end)}
    # {member_name: {native_pos: expected_base}} — what this member has at each discriminating position
    discriminating_positions: dict[str, dict[int, str]]
    # {member_name: {native_pos: (other_members_bases, ...)}} — bases other members have at the same MSA
    # column, excluding gaps and this member's own base. Used to tell MIX from MUT.
    alternative_bases: dict[str, dict[int, tuple[str, ...]]]
    size_class: Literal["point_or_short", "long", "insertion"]


@dataclass
class Group:
    members: tuple[str, ...]          # reference names; size 1 for singletons
    backbone_intervals: dict[str, tuple[tuple[int, int], ...]]  # {ref_name: ((start,end),...)}
    variant_regions: tuple[VariantRegion, ...]


@dataclass
class WindowResult:
    contig: str
    start: int
    end: int
    mean_coverage: float
    status: Literal["intact", "warn", "broken"]
    affected_features: list[str] = field(default_factory=list)


@dataclass
class RegionResult:
    region_index: int
    positions: list[dict]    # [{"pos": int, "A": int, "T": int, "C": int, "G": int, "status": str}]
    region_status: Literal["trusted_clean", "thin_clean", "partial", "mixed", "mutated"]


@dataclass
class MutationHit:
    pos: int        # native coord (0-based)
    expected: str   # reference base at this position
    A: int
    C: int
    G: int
    T: int
    majority_base: str
    majority_fraction: float
    affected_features: list[str] = field(default_factory=list)


@dataclass
class InsertScanResult:
    mode: Literal["quick", "insert", "whole"]
    scan_start: int        # native coord, inclusive
    scan_end: int          # native coord, exclusive
    positions_scanned: int
    mean_coverage: float
    mutations: list[MutationHit] = field(default_factory=list)


@dataclass
class WellResult:
    well_id: str
    total_reads: int
    mapped_reads: int
    winning_group: str | None
    called_member: str | None
    coverage_span: float
    variant_regions: list[RegionResult] = field(default_factory=list)
    backbone_windows: list[WindowResult] = field(default_factory=list)
    insert_scan: InsertScanResult | None = None
    verdict: Literal["GREEN", "YELLOW", "RED_MIX", "RED_MUT",
                     "RED_BROKEN", "RED_UNKNOWN", "NO_DATA"] = "NO_DATA"
    verdict_reason: str = ""
