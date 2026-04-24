from __future__ import annotations

from pathlib import Path

from .types import Feature, Reference


def load_folder(path: Path | str) -> list[Reference]:
    path = Path(path)
    results: list[Reference] = []

    for fpath in sorted(path.iterdir()):
        if fpath.suffix == ".gb":
            results.append(_load_genbank(fpath))
        elif fpath.suffix == ".dna":
            results.append(_load_snapgene(fpath))

    return sorted(results, key=lambda r: r.name)


def _load_genbank(fpath: Path) -> Reference:
    try:
        from Bio import SeqIO
        record = SeqIO.read(str(fpath), "genbank")
        sequence = str(record.seq).upper()

        features: list[Feature] = []
        for feat in record.features:
            if feat.type == "source":
                continue
            label_list = feat.qualifiers.get(
                "label",
                feat.qualifiers.get(
                    "gene",
                    feat.qualifiers.get("product", ["unknown"])
                )
            )
            label = label_list[0]
            start = int(feat.location.start)
            end = int(feat.location.end)
            features.append(Feature(kind=feat.type, label=label, start=start, end=end))

        return Reference(
            name=fpath.stem,
            sequence=sequence,
            features=tuple(features),
            source_path=fpath,
        )
    except Exception as e:
        raise ValueError(f"Failed to load {fpath}: {e}") from e


def _load_snapgene(fpath: Path) -> Reference:
    try:
        import snapgene_reader
        data = snapgene_reader.parse(str(fpath))

        sequence = ""
        if isinstance(data, dict):
            sequence = data.get("seq", data.get("sequence", "")).upper()
        else:
            sequence = str(data).upper()

        features: list[Feature] = []
        try:
            raw_features = data.get("features", []) if isinstance(data, dict) else []
            for feat in raw_features:
                kind = feat.get("type", "misc_feature")
                label = feat.get("label", feat.get("name", "unknown"))
                start = int(feat.get("start", 0))
                end = int(feat.get("end", 0))
                features.append(Feature(kind=kind, label=label, start=start, end=end))
        except Exception:
            features = []

        return Reference(
            name=fpath.stem,
            sequence=sequence,
            features=tuple(features),
            source_path=fpath,
        )
    except Exception as e:
        raise ValueError(f"Failed to load {fpath}: {e}") from e
