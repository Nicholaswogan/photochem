"""Bounded search and reading over the Photochem and model-data repositories."""

from __future__ import annotations

import subprocess
from pathlib import Path

import h5py


ROOT = Path(__file__).resolve().parents[2]
SOURCE_URL = "https://github.com/Nicholaswogan/photochem/blob"
DATA_SOURCE_URL = "https://github.com/Nicholaswogan/photochem_clima_data/blob"
DATA_PREFIX = "photochem_clima_data/"
TEXT_SUFFIXES = {
    ".c", ".cc", ".cmake", ".cpp", ".css", ".f", ".f90", ".h",
    ".html", ".js", ".json", ".md", ".pxd", ".py", ".pyx",
    ".rst", ".toml", ".txt", ".yaml", ".yml",
}
MAX_FILE_BYTES = 1_000_000


class Repository:
    def __init__(self, root: Path = ROOT, data_root: Path | None = None,
                 source_commit: str | None = None, data_commit: str | None = None,
                 release: str | None = None, data_version: str | None = None):
        self.root = root.resolve()
        self.release = release
        self.data_version = data_version
        self.files: dict[str, list[str]] = {}
        self.hdf5_files: dict[str, Path] = {}
        self.sources: dict[str, tuple[str, str, bool]] = {}
        self.commit = self._index(self.root, "", SOURCE_URL, source_commit)
        if data_root is None:
            data_root = self.root.parent / "photochem_clima_data"
        self.data_root = data_root.resolve()
        self.data_commit = None
        if data_commit is not None or (self.data_root / ".git").exists():
            self.data_commit = self._index(self.data_root, DATA_PREFIX, DATA_SOURCE_URL,
                                           data_commit)

    def _index(self, root: Path, prefix: str, source_url: str,
               revision: str | None = None) -> str:
        if revision is None:
            revision = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=root, text=True
            ).strip()
            changed = subprocess.check_output(
                ["git", "diff", "--name-only", "-z", "HEAD"], cwd=root
            )
            changed_paths = {raw.decode("utf-8") for raw in changed.split(b"\0") if raw}
            tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=root)
            paths = [Path(raw.decode("utf-8")) for raw in tracked.split(b"\0") if raw]
        else:
            changed_paths = set()
            paths = [file.relative_to(root) for file in root.rglob("*") if file.is_file()]
        for relative in paths:
            path = relative.as_posix()
            file = (root / relative).resolve()
            if not file.is_relative_to(root) or not file.is_file():
                continue
            name = prefix + path
            if prefix and file.suffix.lower() == ".h5":
                self.hdf5_files[name] = file
                self.sources[name] = (source_url, revision, path in changed_paths)
                continue
            if file.suffix.lower() not in TEXT_SUFFIXES and file.name not in {
                "CMakeLists.txt", "Makefile", "LICENSE",
            }:
                continue
            if file.stat().st_size > MAX_FILE_BYTES:
                continue
            try:
                self.files[name] = file.read_text(encoding="utf-8").splitlines()
                self.sources[name] = (source_url, revision, path in changed_paths)
            except UnicodeDecodeError:
                continue
        return revision

    def source_url(self, path: str, line: int | None = None) -> str:
        source_url, commit, changed = self.sources[path]
        if changed:
            return ""  # An uncommitted line may not exist at the commit URL.
        relative_path = path.removeprefix(DATA_PREFIX) if source_url == DATA_SOURCE_URL else path
        url = f"{source_url}/{commit}/{relative_path}"
        return f"{url}#L{line}" if line is not None else url

    def list_files(self, prefix: str) -> dict:
        if len(prefix) > 100 or prefix.startswith("/") or ".." in Path(prefix).parts:
            return {"error": "Invalid prefix"}
        matches = [path for path in sorted(self.files.keys() | self.hdf5_files.keys())
                   if path.startswith(prefix)]
        return {"files": matches[:100], "total": len(matches)}

    def search_text(self, query: str) -> dict:
        if not query or len(query) > 120:
            return {"error": "Search query must be 1–120 characters"}
        needle = query.casefold()
        matches = []
        for path, lines in sorted(self.files.items()):
            for number, line in enumerate(lines, 1):
                if needle in line.casefold():
                    matches.append({
                        "path": path,
                        "line": number,
                        "text": line[:300],
                        "url": self.source_url(path, number),
                    })
                    if len(matches) == 30:
                        return {"matches": matches, "truncated": True}
        return {"matches": matches, "truncated": False}

    def read_file(self, path: str, start_line: int, end_line: int) -> dict:
        lines = self.files.get(path)
        if lines is None:
            return {"error": "File is not in the tracked text-file allowlist"}
        if not isinstance(start_line, int) or not isinstance(end_line, int):
            return {"error": "Line numbers must be integers"}
        if start_line < 1 or end_line < start_line or end_line - start_line >= 100:
            return {"error": "Request a valid range of at most 100 lines"}
        selected = lines[start_line - 1:end_line]
        return {
            "path": path,
            "start_line": start_line,
            "end_line": start_line + len(selected) - 1,
            "total_lines": len(lines),
            "content": "\n".join(
                f"{number}: {line[:500]}"
                for number, line in enumerate(selected, start_line)
            ),
            "url": self.source_url(path, start_line),
        }

    def inspect_hdf5(self, path: str) -> dict:
        file = self.hdf5_files.get(path)
        if file is None:
            return {"error": "File is not in the tracked HDF5 allowlist"}
        if not file.resolve().is_relative_to(self.data_root) or not file.is_file():
            return {"error": "HDF5 file is outside the data repository"}
        entries = []
        try:
            with h5py.File(file, "r") as data:
                def describe(name, item):
                    if len(entries) >= 50:
                        return True
                    if isinstance(item, h5py.Dataset):
                        entries.append({"name": name, "type": "dataset",
                                        "shape": item.shape, "dtype": str(item.dtype)})
                    else:
                        entries.append({"name": name, "type": "group"})
                    return None

                data.visititems(describe)
        except (OSError, ValueError) as error:
            return {"error": f"Could not inspect HDF5 file: {type(error).__name__}"}
        return {"path": path, "url": self.source_url(path), "entries": entries,
                "truncated": len(entries) == 50}
