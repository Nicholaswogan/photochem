"""Focused checks for the local assistant's repository boundary and tool loop."""

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from fastapi.testclient import TestClient

from docs.ask_ai.release import latest_release_repository, resolve_release
from docs.ask_ai.repository import Repository
from docs.ask_ai.server import (MAX_OUTPUT_TOKENS, OutputTokenLimitError,
                           create_app, run_tool, stream_answer)


class RepositoryTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name) / "photochem"
        self.root.mkdir()
        import subprocess

        subprocess.run(["git", "init", "-q"], cwd=self.root, check=True)
        (self.root / "docs").mkdir()
        (self.root / "docs" / "guide.md").write_text("Photochem guide\nEvoAtmosphere setup\n")
        (self.root / "secret.txt").write_text("untracked secret\n")
        subprocess.run(["git", "add", "docs/guide.md"], cwd=self.root, check=True)
        subprocess.run(
            ["git", "-c", "user.name=Test", "-c", "user.email=test@example.com",
             "commit", "-qm", "Test"], cwd=self.root, check=True
        )
        self.repository = Repository(self.root)

    def tearDown(self):
        self.directory.cleanup()

    def test_only_tracked_text_is_accessible(self):
        self.assertEqual(self.repository.list_files("")["files"], ["docs/guide.md"])
        self.assertEqual(self.repository.search_text("evoatmosphere")["matches"][0]["line"], 2)
        self.assertIn("EvoAtmosphere", self.repository.read_file("docs/guide.md", 1, 2)["content"])
        self.assertIn("error", self.repository.read_file("secret.txt", 1, 2))
        self.assertIn("error", self.repository.read_file("../secret.txt", 1, 2))
        self.assertIn("error", self.repository.read_file("docs/guide.md", 1, 101))

    def test_tool_arguments_cannot_add_capabilities(self):
        self.assertIn("error", run_tool(self.repository, "run_shell", '{"command":"ls"}'))
        self.assertIn("error", run_tool(self.repository, "read_file", "{}"))

    def test_sibling_data_repo_is_searchable_and_hdf5_is_bounded(self):
        import h5py
        import subprocess

        data_root = self.root.parent / "photochem_clima_data"
        data_root.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=data_root, check=True)
        (data_root / "README.md").write_text("Clima reference data\n")
        (data_root / "spectrum.h5").touch()
        with h5py.File(data_root / "spectrum.h5", "w") as data:
            data.create_dataset("wavelength", data=[100.0, 200.0])
        (data_root / "private.txt").write_text("not tracked\n")
        subprocess.run(["git", "add", "README.md", "spectrum.h5"], cwd=data_root, check=True)
        subprocess.run(
            ["git", "-c", "user.name=Test", "-c", "user.email=test@example.com",
             "commit", "-qm", "Data"], cwd=data_root, check=True
        )
        repository = Repository(self.root)
        names = repository.list_files("photochem_clima_data/")["files"]
        self.assertEqual(names, ["photochem_clima_data/README.md",
                                 "photochem_clima_data/spectrum.h5"])
        self.assertIn("Clima", repository.read_file(names[0], 1, 1)["content"])
        self.assertIn("Nicholaswogan/photochem_clima_data/blob",
                      repository.source_url(names[0], 1))
        self.assertEqual(repository.inspect_hdf5(names[1])["entries"][0]["shape"], (2,))
        self.assertIn("error", repository.read_file(names[1], 1, 1))
        self.assertIn("error", repository.inspect_hdf5("photochem_clima_data/private.txt"))


class ToolLoopTests(unittest.TestCase):
    def test_token_exhaustion_is_reported(self):
        client = Mock()
        client.responses.create.return_value = iter([
            SimpleNamespace(
                type="response.incomplete",
                response=SimpleNamespace(
                    incomplete_details=SimpleNamespace(reason="max_output_tokens")
                ),
            ),
        ])

        with self.assertRaises(OutputTokenLimitError):
            list(stream_answer(client, Mock(), "Explain the model", []))
        request = client.responses.create.call_args.kwargs
        self.assertEqual(request["reasoning"]["effort"], "medium")
        self.assertEqual(request["max_output_tokens"], MAX_OUTPUT_TOKENS)

    def test_function_call_is_returned_to_model_and_answer_streams(self):
        repository = Mock()
        repository.search_text.return_value = {"matches": [{"path": "docs/guide.md"}]}
        first = SimpleNamespace(
            id="response-1",
            output=[SimpleNamespace(type="function_call", name="search_text",
                                    arguments=json.dumps({"query": "guide"}), call_id="call-1")],
            output_text="",
        )
        second = SimpleNamespace(id="response-2", output=[], output_text="Found the guide.")
        client = Mock()
        client.responses.create.side_effect = [
            iter([SimpleNamespace(type="response.completed", response=first)]),
            iter([
                SimpleNamespace(type="response.output_text.delta", delta="Found "),
                SimpleNamespace(type="response.output_text.delta", delta="the guide."),
                SimpleNamespace(type="response.completed", response=second),
            ]),
        ]

        events = list(stream_answer(client, repository, "Where is the guide?", []))

        self.assertEqual(events[-1], {"type": "done", "answer": "Found the guide."})
        self.assertEqual([event["text"] for event in events if event["type"] == "delta"],
                         ["Found ", "the guide."])
        repository.search_text.assert_called_once_with("guide")
        followup = client.responses.create.call_args_list[1].kwargs
        self.assertEqual(followup["previous_response_id"], "response-1")
        self.assertEqual(followup["input"][0]["call_id"], "call-1")

    def test_http_endpoint_returns_ndjson(self):
        repository = Mock(commit="0123456789abcdef", data_commit=None,
                          files=["docs/guide.md"], hdf5_files={},
                          release="v0.9.0", data_version="0.3.2")
        response = SimpleNamespace(id="r1", output=[], output_text="Hello.")
        client = Mock()
        client.responses.create.return_value = iter([
            SimpleNamespace(type="response.output_text.delta", delta="Hello."),
            SimpleNamespace(type="response.completed", response=response),
        ])
        web = TestClient(create_app(repository, client))
        self.assertEqual(web.get("/health").json()["protocol"], "ndjson-v1")
        self.assertEqual(web.get("/health").json()["release"], "v0.9.0")
        result = web.post("/chat", json={"message": "Hello?"})
        self.assertEqual(result.status_code, 200)
        self.assertEqual(result.headers["content-type"], "application/x-ndjson")
        self.assertEqual([json.loads(line)["type"] for line in result.text.splitlines()],
                         ["delta", "done"])
        self.assertEqual(web.post("/chat", json={"message": " "}).status_code, 422)


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name) / "photochem"
        self.data_root = Path(self.directory.name) / "photochem_clima_data"
        for repo in (self.root, self.data_root):
            repo.mkdir()
            subprocess.run(["git", "init", "-q"], cwd=repo, check=True)

        (self.data_root / "pyproject.toml").write_text('[project]\nversion = "0.3.2"\n')
        (self.data_root / "data.txt").write_text("released data\n")
        import h5py
        with h5py.File(self.data_root / "spectrum.h5", "w") as data:
            data.create_dataset("wavelength", data=[100.0, 200.0])
        self.data_sha = self.commit(self.data_root)
        subprocess.run(["git", "tag", "v0.3.2"], cwd=self.data_root, check=True)
        (self.data_root / "data.txt").write_text("unreleased data\n")
        with h5py.File(self.data_root / "spectrum.h5", "w") as data:
            data.create_dataset("newer", data=[1.0])
        self.commit(self.data_root)

        (self.root / "CMakeLists.txt").write_text(
            'project(Photochem LANGUAGES Fortran C VERSION "0.9.0")\n'
            'set(PHOTOCHEM_CLIMA_DATA_VERSION "0.3.2")\n'
        )
        (self.root / "guide.md").write_text("released Photochem\n")
        self.source_sha = self.commit(self.root)
        subprocess.run(["git", "tag", "v0.9.0"], cwd=self.root, check=True)
        (self.root / "guide.md").write_text("unreleased Photochem\n")
        self.commit(self.root)

    def tearDown(self):
        self.directory.cleanup()

    @staticmethod
    def commit(repo):
        subprocess.run(["git", "add", "-A"], cwd=repo, check=True)
        subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.com",
                        "commit", "-qm", "Test"], cwd=repo, check=True)
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo,
                                       text=True).strip()

    def test_indexes_matching_release_commits_only(self):
        self.assertEqual(resolve_release(self.root, self.data_root),
                         ("v0.9.0", self.source_sha, "0.3.2", self.data_sha))
        with latest_release_repository(self.root, self.data_root,
                                       Path(self.directory.name)) as repository:
            self.assertEqual(repository.release, "v0.9.0")
            self.assertEqual(repository.data_version, "0.3.2")
            self.assertIn("released Photochem", repository.read_file("guide.md", 1, 1)["content"])
            self.assertIn("released data", repository.read_file(
                "photochem_clima_data/data.txt", 1, 1)["content"])
            self.assertNotIn("unreleased", repository.read_file("guide.md", 1, 1)["content"])
            self.assertIn(self.source_sha, repository.source_url("guide.md", 1))
            self.assertIn(self.data_sha, repository.source_url(
                "photochem_clima_data/data.txt", 1))
            entries = repository.inspect_hdf5("photochem_clima_data/spectrum.h5")["entries"]
            self.assertEqual([entry["name"] for entry in entries], ["wavelength"])

    def test_selects_newest_version_tag(self):
        (self.root / "CMakeLists.txt").write_text(
            'project(Photochem LANGUAGES Fortran C VERSION "0.10.0")\n'
            'set(PHOTOCHEM_CLIMA_DATA_VERSION "0.3.2")\n'
        )
        newer_sha = self.commit(self.root)
        subprocess.run(["git", "tag", "v0.10.0"], cwd=self.root, check=True)
        self.assertEqual(resolve_release(self.root, self.data_root)[:2],
                         ("v0.10.0", newer_sha))

    def test_missing_data_tag_fails_instead_of_using_head(self):
        subprocess.run(["git", "tag", "-d", "v0.3.2"], cwd=self.data_root,
                       check=True, capture_output=True)
        with self.assertRaisesRegex(RuntimeError, "Missing photochem_clima_data tag"):
            resolve_release(self.root, self.data_root)

if __name__ == "__main__":
    unittest.main()
