from __future__ import annotations

import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from Bio import SeqIO

import bm_pirna.download_protein_sequences as downloader


class FakeResponse:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *args: object) -> None:
        return None

    def read(self) -> bytes:
        return self.payload


class DownloadProteinSequencesTests(unittest.TestCase):
    def test_single_column_bom_crlf_and_mixed_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ids.txt"
            path.write_bytes(b"\xef\xbb\xbfNP_001234.1\r\nFBpp0000123\r\nFBgn0000456\r\n")
            parsed = downloader.parse_input_file(path)
        self.assertEqual(parsed.id_column, "column_1")
        self.assertEqual(len(parsed.entries), 3)
        self.assertEqual(
            [entry.detected_type for entry in parsed.entries],
            [
                "ncbi_protein_accession",
                "flybase_protein_id",
                "flybase_gene_id",
            ],
        )

    def test_single_column_symbol_is_not_treated_as_header(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "symbols.txt"
            path.write_text("CREB\nRack1\n", encoding="utf-8")
            parsed = downloader.parse_input_file(path)
        self.assertEqual(len(parsed.entries), 2)
        self.assertEqual(
            [entry.request_kind for entry in parsed.entries],
            ["ncbi_symbol", "ncbi_symbol"],
        )

    def test_header_alias_and_column_priority(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.tsv"
            path.write_text(
                "gene_name\tgene_id\tprotein_id\nCrebB\t692871\tNP_001234.1\n",
                encoding="utf-8",
            )
            parsed = downloader.parse_input_file(path)
        self.assertEqual(parsed.id_column, "protein_id")
        self.assertEqual(parsed.entries[0].request_kind, "ncbi_protein")

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genes.tsv"
            path.write_text("gene_name\tgene_id\nCrebB\t692871\n", encoding="utf-8")
            parsed = downloader.parse_input_file(path)
        self.assertEqual(parsed.id_column, "gene_id")
        self.assertEqual(parsed.entries[0].request_id, "692871")

    def test_csv_and_explicit_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.csv"
            path.write_text(
                "sample,identifier\nA,FBgn0000001\nB,FBgn0000002\n",
                encoding="utf-8",
            )
            parsed = downloader.parse_input_file(path, id_column="identifier")
        self.assertEqual(parsed.id_column, "identifier")
        self.assertEqual(
            [entry.request_id for entry in parsed.entries],
            [
                "FBgn0000001",
                "FBgn0000002",
            ],
        )

    def test_classification_respects_explicit_type_and_source(self) -> None:
        self.assertEqual(
            downloader.classify_id("123", id_type="auto"),
            ("ncbi_gene_id", "ncbi", "123"),
        )
        self.assertEqual(
            downloader.classify_id("NP_001234.1", id_type="gene")[0],
            "unrecognized",
        )
        self.assertEqual(
            downloader.classify_id("FBpp0000001", id_type="gene")[0],
            "unrecognized",
        )
        self.assertEqual(
            downloader.classify_id("FBgn0000001", id_type="protein")[0],
            "unsupported_gene_protein",
        )
        self.assertEqual(
            downloader.classify_id("123", source="flybase")[0],
            "source_mismatch",
        )

    def test_gene_map_resolution_and_unmapped_loc(self) -> None:
        entries = [
            downloader.InputEntry(
                input_file="x",
                row_number=1,
                id_column="gene_name",
                raw_id="CREB",
                normalized_id="CREB",
                detected_type="ncbi_gene_symbol",
                source="ncbi",
                request_kind="ncbi_symbol",
                request_id="CREB",
            ),
            downloader.InputEntry(
                input_file="x",
                row_number=2,
                id_column="gene_name",
                raw_id="LOC123",
                normalized_id="LOC123",
                detected_type="ncbi_gene_symbol",
                source="ncbi",
                request_kind="ncbi_symbol",
                request_id="LOC123",
            ),
        ]
        downloader._resolve_symbol_entries(
            entries,
            gene_map={"creb": {"692871"}},
            taxid=None,
        )
        self.assertEqual(entries[0].request_kind, "ncbi_gene")
        self.assertEqual(entries[0].request_id, "692871")
        self.assertEqual(entries[1].status, "ambiguous")
        self.assertIn("taxid", entries[1].message)

    def test_symbol_resolution_uses_exact_taxon_results(self) -> None:
        stdout = (
            '{"gene_id":"692871","symbol":"CrebB","synonyms":["CREB"]}\n'
            '{"gene_id":"999999","symbol":"Other","synonyms":[]}\n'
        )
        completed = downloader.subprocess.CompletedProcess(
            ["datasets"], 0, stdout=stdout, stderr=""
        )
        with patch.object(downloader, "_run_command", return_value=completed) as run:
            result = downloader.resolve_symbols_with_datasets(["CREB"], taxid="7091")
        self.assertEqual(result, {"CREB": ("692871", "")})
        self.assertEqual(run.call_count, 1)
        command = run.call_args.args[0]
        self.assertIn("--taxon", command)
        self.assertIn("7091", command)

    def test_choose_longest_and_deterministic_tie(self) -> None:
        proteins = [
            downloader.ProteinSequence("B", "MKTK"),
            downloader.ProteinSequence("A", "AAAA"),
            downloader.ProteinSequence("C", "MKT"),
        ]
        selected = downloader._choose_proteins(proteins, "longest")
        self.assertEqual([protein.accession for protein in selected], ["A"])
        self.assertEqual(
            [protein.accession for protein in downloader._choose_proteins(proteins, "all")],
            ["A", "B", "C"],
        )

    def test_nested_flybase_payload_and_sequence_cleanup(self) -> None:
        payload = {
            "resultset": {
                "data_version": "FB2026_01",
                "result": [
                    {
                        "id": "FBpp0000002",
                        "sequence": "M K\nT",
                        "description": "short",
                    },
                    {
                        "id": "FBpp0000001",
                        "sequence": "MKTK",
                        "description": "long",
                    },
                ],
                "not_found": [],
            }
        }
        proteins, version, message = downloader._proteins_from_flybase_payload(
            payload, source_gene_id="FBgn0000001"
        )
        self.assertEqual(version, "FB2026_01")
        self.assertEqual(message, "")
        self.assertEqual([len(protein.sequence) for protein in proteins], [3, 4])
        self.assertTrue(all(protein.gene_id == "FBgn0000001" for protein in proteins))

    def test_flybase_retry_on_rate_limit(self) -> None:
        payload = json.dumps(
            {
                "resultset": {
                    "data_version": "FB2026_01",
                    "result": [{"id": "FBpp1", "sequence": "MKT"}],
                }
            }
        ).encode()
        calls = {"count": 0}

        def fake_urlopen(request: object, timeout: float) -> FakeResponse:
            calls["count"] += 1
            if calls["count"] == 1:
                raise downloader.HTTPError(
                    "https://api.flybase.org", 429, "rate limited", {}, io.BytesIO()
                )
            return FakeResponse(payload)

        with (
            patch.object(downloader, "urlopen", side_effect=fake_urlopen),
            patch.object(downloader.time, "sleep"),
        ):
            result = downloader._request_flybase(
                "FBgn0000001",
                subtype="FBpp",
                timeout=1,
                retries=2,
                limiter=downloader.RateLimiter(0),
            )
        self.assertEqual(calls["count"], 2)
        self.assertEqual(result["resultset"]["result"][0]["id"], "FBpp1")

    def test_ncbi_batch_recursion_preserves_accession_results(self) -> None:
        def fake_run(command: list[str]) -> downloader.subprocess.CompletedProcess[str]:
            input_file = Path(command[command.index("--inputfile") + 1])
            package = Path(command[command.index("--filename") + 1])
            accessions = input_file.read_text(encoding="utf-8").split()
            if len(accessions) > 1:
                return downloader.subprocess.CompletedProcess(
                    command, 1, stdout="", stderr="batch failed"
                )
            accession = accessions[0]
            with zipfile.ZipFile(package, "w") as archive:
                archive.writestr(
                    "ncbi_dataset/data/protein.faa",
                    f">{accession} model [GeneID=123]\nMKT\n",
                )
            return downloader.subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        with patch.object(downloader, "_run_command", side_effect=fake_run):
            proteins, _, _, errors = downloader._run_ncbi_batches(
                ["NP_000001.1", "NP_000002.1"],
                mode="accession",
                batch_size=2,
                retries=1,
            )
        self.assertEqual(errors, {})
        self.assertEqual(set(proteins), {"np_000001.1", "np_000002.1"})

    def test_ncbi_package_parser_uses_expected_members(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            package = Path(tmp) / "package.zip"
            with zipfile.ZipFile(package, "w") as archive:
                archive.writestr("../data/protein.faa", ">BAD\nBAD\n")
                archive.writestr(
                    "ncbi_dataset/data/data_report.jsonl",
                    json.dumps(
                        {
                            "geneId": "123",
                            "annotations": [{"annotationName": "RS_2026"}],
                        }
                    )
                    + "\n",
                )
                archive.writestr(
                    "ncbi_dataset/data/protein.faa",
                    ">NP_000001.1 protein [GeneID=123]\nMKT\n",
                )
            proteins, gene_ids, versions = downloader._parse_ncbi_package(package)
        self.assertEqual(gene_ids, {"123"})
        self.assertEqual(versions, {"123": "RS_2026"})
        self.assertEqual(proteins["123"][0].accession, "NP_000001.1")

    def test_ncbi_efetch_fallback_parses_direct_accessions(self) -> None:
        payload = b">XP_004929423.1 model protein [GeneID=123]\nMKT\n"
        calls = {"count": 0}

        def fake_urlopen(request: object, timeout: float) -> FakeResponse:
            calls["count"] += 1
            self.assertIn("XP_004929423.1", request.full_url)
            return FakeResponse(payload)

        with patch.object(downloader, "urlopen", side_effect=fake_urlopen):
            proteins, errors = downloader._fetch_ncbi_proteins_efetch(
                ["XP_004929423.1"], timeout=1, retries=1, limiter=downloader.RateLimiter(0)
            )
        self.assertEqual(errors, {})
        self.assertEqual(calls["count"], 1)
        self.assertEqual(proteins["xp_004929423.1"][0].gene_id, "123")

    def test_run_download_writes_outputs_and_preserves_duplicate_rows(self) -> None:
        def fake_ncbi(
            request_ids: list[str],
            *,
            mode: str,
            batch_size: int,
            retries: int,
        ) -> tuple[
            dict[str, list[downloader.ProteinSequence]],
            set[str],
            dict[str, str],
            dict[str, str],
        ]:
            if mode == "gene":
                return (
                    {
                        "123": [
                            downloader.ProteinSequence(
                                "NP_GENE", "MKT", source="NCBI", gene_id="123"
                            )
                        ]
                    },
                    {"123"},
                    {"123": "RS_2026"},
                    {},
                )
            return (
                {
                    "np_000002.1": [
                        downloader.ProteinSequence(
                            "NP_000002.1", "MK", source="NCBI", gene_id="456"
                        )
                    ]
                },
                {"456"},
                {},
                {},
            )

        def fake_flybase(
            entries: object,
            *,
            timeout: float,
            retries: int,
            limiter: object = None,
        ) -> dict[tuple[str, str], downloader.Resolution]:
            return {
                ("flybase_gene", "FBgn0000001"): downloader.Resolution(
                    status="downloaded",
                    proteins=[
                        downloader.ProteinSequence(
                            "FBpp_LONG", "MKTK", source="FlyBase", gene_id="FBgn0000001"
                        )
                    ],
                    resolved_gene_id="FBgn0000001",
                    data_version="FB2026_01",
                )
            }

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gene_file = root / "genes.tsv"
            gene_file.write_text("gene_id\n123\n123\n", encoding="utf-8")
            protein_file = root / "proteins.txt"
            protein_file.write_text("NP_000002.1\n", encoding="utf-8")
            fly_file = root / "fly.txt"
            fly_file.write_text("FBgn0000001\n", encoding="utf-8")
            output = root / "out"
            with (
                patch.object(downloader, "_run_ncbi_batches", side_effect=fake_ncbi),
                patch.object(downloader, "fetch_flybase_resolutions", side_effect=fake_flybase),
                patch.object(downloader, "datasets_version", return_value="datasets 18.36.0"),
            ):
                metadata = downloader.run_download(
                    [gene_file, protein_file, fly_file],
                    output_dir=output,
                    force=False,
                )
            with (output / "all.proteins.faa").open(encoding="utf-8") as fasta_handle:
                records = list(SeqIO.parse(fasta_handle, "fasta"))
            manifest = (
                (output / "protein_download_manifest.tsv").read_text(encoding="utf-8").splitlines()
            )
            stored = json.loads((output / "download_metadata.json").read_text())
        self.assertEqual(metadata["counts"]["input_rows"], 4)
        self.assertEqual(metadata["counts"]["unique_proteins"], 3)
        self.assertEqual(len(records), 3)
        self.assertEqual(len(manifest), 5)
        self.assertEqual(stored["ncbi_datasets_version"], "datasets 18.36.0")
        self.assertEqual(stored["flybase_data_versions"], ["FB2026_01"])

    def test_allow_partial_writes_manifest_for_invalid_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "invalid.txt"
            path.write_text("not-an-id\n", encoding="utf-8")
            output = root / "out"
            metadata = downloader.run_download([path], output_dir=output, allow_partial=True)
            manifest_exists = (output / "protein_download_manifest.tsv").exists()
        self.assertEqual(metadata["counts"]["fatal_rows"], 1)
        self.assertTrue(manifest_exists)

    def test_empty_id_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "empty.txt"
            path.write_text("NA\n", encoding="utf-8")
            output = root / "out"
            metadata = downloader.run_download([path], output_dir=output, allow_partial=True)
        self.assertEqual(metadata["counts"]["fatal_rows"], 1)


if __name__ == "__main__":
    unittest.main()
