"""Textual TUI for pickme.

Main screen: list of runs found under `input/` with status chips.
Run screen: live per-well progress while `run_pipeline` executes.
"""

from __future__ import annotations

import threading
from pathlib import Path
from typing import Optional

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import Screen
from textual.widgets import DataTable, Footer, Header, Static

from pickme.runs import Run, scan_input_dir
from pickme.runner import run_pipeline


DEFAULT_INPUT = Path("input")
DEFAULT_OUTPUT = Path("output")


_STATUS_STYLE = {
    "new": "bold #007AFF",
    "done": "bold #1a7a3a",
    "partial": "bold #b26a00",
    "invalid": "bold #c0362c",
}


def _status_markup(run: Run) -> str:
    style = _STATUS_STYLE.get(run.status, "")
    label = f"[{run.status}]"
    if style:
        return f"[{style}]{label}[/]"
    return label


class RunScreen(Screen):
    """Live progress screen for a single run."""

    BINDINGS = [
        Binding("q", "app.pop_screen", "Back"),
        Binding("escape", "app.pop_screen", "Back"),
    ]

    def __init__(self, run: Run, scan_mode: str = "quick") -> None:
        super().__init__()
        self._run = run
        self._scan_mode = scan_mode
        self._table: Optional[DataTable] = None
        self._header: Optional[Static] = None
        self._footer: Optional[Static] = None
        self._row_keys: dict[str, str] = {}
        self._worker: Optional[threading.Thread] = None

    def compose(self) -> ComposeResult:
        yield Header()
        self._header = Static(
            f"Running pickme on [b]{self._run.name}[/b] (mode: {self._scan_mode})",
            id="run-header",
        )
        yield self._header
        table: DataTable = DataTable(id="well-table")
        table.add_columns("well", "stage", "detail")
        self._table = table
        yield table
        self._footer = Static("Starting...", id="run-footer")
        yield self._footer
        yield Footer()

    def on_mount(self) -> None:
        self._worker = threading.Thread(target=self._work, daemon=True)
        self._worker.start()

    # --- worker side ---

    def _work(self) -> None:
        try:
            run_pipeline(
                self._run,
                scan_mode=self._scan_mode,
                on_progress=self._on_progress_threadsafe,
            )
        except Exception as exc:  # noqa: BLE001
            self.app.call_from_thread(self._render_error, exc)

    def _on_progress_threadsafe(
        self, well_id: str, stage: str, payload: dict
    ) -> None:
        self.app.call_from_thread(self._on_progress, well_id, stage, payload)

    # --- UI thread side ---

    def _on_progress(self, well_id: str, stage: str, payload: dict) -> None:
        assert self._table is not None
        assert self._footer is not None

        if stage == "start":
            self._footer.update(
                f"Aligning {payload['n_wells']} wells..."
            )
            return

        if stage == "summary":
            verdicts = payload.get("verdicts", {})
            parts = [f"{v}: {c}" for v, c in sorted(verdicts.items())]
            self._footer.update(
                "Done. " + "  ".join(parts)
                + f"   [q] back"
            )
            return

        if not well_id:
            return

        if stage == "align":
            detail = f"{payload.get('mapped_reads', 0)}/{payload.get('total_reads', 0)} mapped"
        elif stage == "identify":
            verdict = payload.get("verdict", "?")
            called = payload.get("called_member", "unknown")
            flagged = payload.get("flagged", 0)
            detail = f"{verdict} · {called}"
            if flagged:
                detail += f" · {flagged} flagged"
        else:
            detail = str(payload)

        if well_id in self._row_keys:
            row_key = self._row_keys[well_id]
            self._table.update_cell(row_key, "stage", stage)
            self._table.update_cell(row_key, "detail", detail)
        else:
            row_key = self._table.add_row(well_id, stage, detail, key=well_id)
            self._row_keys[well_id] = row_key

    def _render_error(self, exc: Exception) -> None:
        assert self._footer is not None
        self._footer.update(
            f"[bold red]Error:[/bold red] {exc}   [q] back"
        )


class PickmeApp(App):
    """Main TUI app: list runs, run one on enter."""

    CSS = """
    Screen {
        background: #ffffff;
        color: #1c1c1e;
    }
    Header {
        background: #f5f5f7;
        color: #1c1c1e;
    }
    Footer {
        background: #f5f5f7;
        color: #1c1c1e;
    }
    #run-header, #run-footer, #input-label {
        padding: 0 1;
        background: #ffffff;
        color: #1c1c1e;
    }
    DataTable {
        height: 1fr;
        background: #ffffff;
        color: #1c1c1e;
    }
    DataTable > .datatable--header {
        background: #f5f5f7;
        color: #6e6e73;
    }
    DataTable > .datatable--cursor {
        background: #e8e8ed;
        color: #1c1c1e;
    }
    DataTable > .datatable--hover {
        background: #f5f5f7;
    }
    """

    TITLE = "pickme"
    SUB_TITLE = "tn5 verification"

    BINDINGS = [
        Binding("enter", "run_selected", "Run"),
        Binding("r", "refresh", "Refresh"),
        Binding("q", "quit", "Quit"),
    ]

    def __init__(
        self,
        input_dir: Path = DEFAULT_INPUT,
        output_dir: Path = DEFAULT_OUTPUT,
        scan_mode: str = "quick",
    ) -> None:
        super().__init__()
        try:
            self.theme = "textual-light"
        except Exception:
            try:
                self.dark = False
            except Exception:
                pass
        self._input_dir = Path(input_dir).resolve()
        self._output_dir = Path(output_dir).resolve()
        self._scan_mode = scan_mode
        self._runs: list[Run] = []
        self._table: Optional[DataTable] = None

    def on_mount(self) -> None:
        self._refresh_table()

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static(
            f"Runs in [b]{self._input_dir}[/b]",
            id="input-label",
        )
        table: DataTable = DataTable(id="runs-table", cursor_type="row")
        table.add_columns("name", "status", "summary")
        self._table = table
        yield table
        yield Footer()

    def _refresh_table(self) -> None:
        assert self._table is not None
        self._table.clear()
        self._runs = scan_input_dir(self._input_dir, self._output_dir)
        if not self._runs:
            self._table.add_row(
                "(no runs)",
                "",
                f"drop a folder into {self._input_dir}",
            )
            return
        for r in self._runs:
            self._table.add_row(
                r.name,
                _status_markup(r),
                r.summary,
            )

    def action_refresh(self) -> None:
        self._refresh_table()

    def action_run_selected(self) -> None:
        assert self._table is not None
        if not self._runs:
            return
        row_index = self._table.cursor_row
        if row_index is None or row_index < 0 or row_index >= len(self._runs):
            return
        run = self._runs[row_index]
        if run.status == "invalid":
            self.bell()
            return
        self.push_screen(RunScreen(run, scan_mode=self._scan_mode))


def main() -> None:  # thin wrapper for console_scripts, though __main__ is canonical
    PickmeApp().run()


if __name__ == "__main__":
    main()
