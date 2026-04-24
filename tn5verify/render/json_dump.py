from __future__ import annotations

import json
import dataclasses
from pathlib import Path


def _default(obj):
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return dataclasses.asdict(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Not serialisable: {type(obj)}")


def write(results, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump([dataclasses.asdict(r) for r in results], f, indent=2, default=_default)


def load(path: Path) -> list[dict]:
    path = Path(path)
    with open(path) as f:
        return json.load(f)
