"""
features.py — M6: annotate affected features on broken/warn backbone windows.
"""
from __future__ import annotations

from .types import Reference, WellResult

_PAD = 300  # bases prepended to sequences >= 500 bp


def annotate_gaps(result: WellResult, refs: list[Reference]) -> None:
    """
    Mutate result in-place.  For each WindowResult in result.backbone_windows
    whose status is "warn" or "broken", find which annotated features overlap
    the window (in native coords) and populate window.affected_features.
    """
    if not result.backbone_windows:
        return
    if result.called_member is None:
        return

    # Find the reference for the called member
    ref: Reference | None = None
    for r in refs:
        if r.name == result.called_member:
            ref = r
            break

    if ref is None:
        return

    seq_len = len(ref.sequence)
    is_padded = seq_len >= 500

    for window in result.backbone_windows:
        if window.status not in ("warn", "broken"):
            continue

        # Convert padded coords back to native coords
        if is_padded:
            native_start = max(0, window.start - _PAD)
            native_end = min(seq_len, window.end - _PAD)
        else:
            native_start = max(0, window.start)
            native_end = min(seq_len, window.end)

        # Find overlapping features
        for feature in ref.features:
            if feature.start < native_end and feature.end > native_start:
                window.affected_features.append(
                    f"{feature.label} ({feature.kind})"
                )

        # Fallback if no features found
        if not window.affected_features:
            window.affected_features.append(f"bp {native_start}\u2013{native_end}")
