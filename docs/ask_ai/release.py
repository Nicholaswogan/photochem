"""Resolve and materialize the latest supported Photochem release pair."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from contextlib import contextmanager
from pathlib import Path, PurePosixPath

from .repository import ROOT, Repository


RELEASE_TAG = re.compile(r"^v(\d+)\.(\d+)\.(\d+)$")
PHOTOCHEM_VERSION = re.compile(
    r'^project\(Photochem\s+LANGUAGES\s+Fortran\s+C\s+VERSION\s+"([^"]+)"\)', re.M
)
DATA_VERSION = re.compile(r'^set\(PHOTOCHEM_CLIMA_DATA_VERSION\s+"([^"]+)"\)', re.M)


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def resolve_release(root: Path = ROOT, data_root: Path | None = None) -> tuple[str, str, str, str]:
    """Return Photochem tag/SHA and matching data version/SHA, or fail closed."""
    root = root.resolve()
    data_root = (data_root or root.parent / "photochem_clima_data").resolve()
    if not (data_root / ".git").exists():
        raise RuntimeError(f"Missing photochem_clima_data Git checkout: {data_root}")
    versions = [(tuple(map(int, match.groups())), tag)
                for tag in git(root, "tag", "--list").splitlines()
                if (match := RELEASE_TAG.fullmatch(tag))]
    if not versions:
        raise RuntimeError("No stable Photochem release tags are available")
    _, tag = max(versions)
    source_sha = git(root, "rev-parse", f"refs/tags/{tag}^{{commit}}")
    cmake = git(root, "show", f"{source_sha}:CMakeLists.txt")
    project = PHOTOCHEM_VERSION.search(cmake)
    data = DATA_VERSION.search(cmake)
    if not project or project.group(1) != tag.removeprefix("v") or not data:
        raise RuntimeError(f"Release {tag} has missing or inconsistent CMake version pins")
    data_version = data.group(1)
    data_tag = f"refs/tags/v{data_version}"
    if subprocess.run(["git", "rev-parse", "--verify", "-q", f"{data_tag}^{{commit}}"],
                      cwd=data_root, capture_output=True).returncode != 0:
        raise RuntimeError(f"Missing photochem_clima_data tag v{data_version}; "
                           "run `git fetch --tags` in that checkout")
    try:
        data_sha = git(data_root, "rev-parse", f"{data_tag}^{{commit}}")
        manifest = git(data_root, "show", f"{data_sha}:pyproject.toml")
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"Data revision for {data_version} is unavailable locally") from error
    actual = re.search(r'^version\s*=\s*"([^"]+)"', manifest, re.M)
    if not actual or actual.group(1) != data_version:
        raise RuntimeError(f"Data revision does not declare version {data_version}")
    return tag, source_sha, data_version, data_sha


def extract_commit(repo: Path, sha: str, destination: Path) -> None:
    destination.mkdir()
    archive = destination.parent / f"{destination.name}.tar"
    try:
        subprocess.run(["git", "archive", "--format=tar", "--output", str(archive), sha],
                       cwd=repo, check=True, capture_output=True)
        with tarfile.open(archive) as tar:
            for member in tar:
                relative = PurePosixPath(member.name)
                if relative.is_absolute() or ".." in relative.parts:
                    raise RuntimeError("Invalid path in release archive")
                target = destination.joinpath(*relative.parts)
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with tar.extractfile(member) as source, target.open("wb") as output:
                        shutil.copyfileobj(source, output)
    finally:
        archive.unlink(missing_ok=True)


@contextmanager
def latest_release_repository(root: Path = ROOT, data_root: Path | None = None,
                              snapshot_parent: Path | None = None):
    root = root.resolve()
    data_root = (data_root or root.parent / "photochem_clima_data").resolve()
    tag, source_sha, data_version, data_sha = resolve_release(root, data_root)
    if snapshot_parent is None:
        snapshot_parent = Path(os.environ.get(
            "ASK_AI_SNAPSHOT_DIR", Path.home() / ".cache" / "photochem-ask-ai"
        ))
    snapshot_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="release-", dir=snapshot_parent) as directory:
        snapshot = Path(directory)
        source_snapshot = snapshot / "photochem"
        data_snapshot = snapshot / "photochem_clima_data"
        extract_commit(root, source_sha, source_snapshot)
        extract_commit(data_root, data_sha, data_snapshot)
        yield Repository(source_snapshot, data_snapshot, source_sha, data_sha,
                         tag, data_version)
