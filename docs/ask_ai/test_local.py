"""Focused checks for the local assistant's repository boundary and tool loop."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from fastapi.testclient import TestClient

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
                          files=["docs/guide.md"], hdf5_files={})
        response = SimpleNamespace(id="r1", output=[], output_text="Hello.")
        client = Mock()
        client.responses.create.return_value = iter([
            SimpleNamespace(type="response.output_text.delta", delta="Hello."),
            SimpleNamespace(type="response.completed", response=response),
        ])
        web = TestClient(create_app(repository, client))
        self.assertEqual(web.get("/health").json()["protocol"], "ndjson-v1")
        result = web.post("/chat", json={"message": "Hello?"})
        self.assertEqual(result.status_code, 200)
        self.assertEqual(result.headers["content-type"], "application/x-ndjson")
        self.assertEqual([json.loads(line)["type"] for line in result.text.splitlines()],
                         ["delta", "done"])
        self.assertEqual(web.post("/chat", json={"message": " "}).status_code, 422)

if __name__ == "__main__":
    unittest.main()
