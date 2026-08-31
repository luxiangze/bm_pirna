from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from bm_pirna import find_homolog_pairs as homolog


def hit_line(query: str, subject: str, **values: object) -> str:
    defaults: dict[str, object] = {
        "identity": 80.0,
        "length": 100,
        "qlen": 100,
        "slen": 100,
        "qcov": 100.0,
        "scov": 100.0,
        "evalue": 1e-20,
        "bitscore": 200.0,
    }
    defaults.update(values)
    fields = [
        query,
        subject,
        str(defaults["identity"]),
        str(defaults["length"]),
        str(defaults["qlen"]),
        str(defaults["slen"]),
        str(defaults["qcov"]),
        str(defaults["scov"]),
        str(defaults["evalue"]),
        str(defaults["bitscore"]),
    ]
    return "\t".join(fields)


class HomologPairTests(unittest.TestCase):
    def test_example_headers_map_every_protein(self) -> None:
        paths = [
            Path(
                "data/external/homology/protein_sequences/bm_AS_genes_1_control_vs_sugp1.proteins.faa"
            ),
            Path(
                "data/external/homology/protein_sequences/dm_AS_genes_control_vs_sugp1.proteins.faa"
            ),
        ]
        for path in paths:
            records = homolog._read_fasta(path)
            mapping = homolog._mapping_from_headers(records)
            self.assertEqual(len(records), len(mapping))
            self.assertTrue(all(value.gene_id for value in mapping.values()))

    def test_mapping_file_supports_header_and_labels(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "map.tsv"
            path.write_text(
                "protein_id\tgene_id\tlabel\nP1\t00123\tLOC123\nP2\tFBgn0000001\tgene-x\n",
                encoding="utf-8",
            )
            mapping = homolog._read_mapping_file(path)
        self.assertEqual(mapping["P1"], homolog.ProteinMapping("123", "LOC123"))
        self.assertEqual(mapping["P2"].gene_id, "FBgn0000001")

    def test_manifest_mapping_filters_by_fasta_stem(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "sample.proteins.faa"
            fasta.write_text(">P1\nMKT\n", encoding="utf-8")
            manifest = root / "protein_download_manifest.tsv"
            manifest.write_text(
                "input_file\tprotein_accession\tresolved_gene_id\traw_id\tstatus\n"
                f"{root / 'sample.txt'}\tP1\t123\tLOC123\tdownloaded\n"
                f"{root / 'other.txt'}\tP1\t999\tOTHER\tdownloaded\n",
                encoding="utf-8",
            )
            mapping = homolog._manifest_mapping(manifest, fasta)
        self.assertEqual(mapping, {"P1": homolog.ProteinMapping("123", "LOC123")})

    def test_parse_hits_and_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.tsv"
            path.write_text(
                hit_line("P1", "Q1", identity=40, qcov=70, scov=70)
                + "\n"
                + hit_line("P1", "Q2", identity=39.9)
                + "\n",
                encoding="utf-8",
            )
            hits = homolog._parse_hits(path)
        self.assertEqual(len(hits), 2)
        filtered = homolog._filter_hits(
            hits,
            min_identity=40,
            min_query_cover=70,
            min_subject_cover=70,
            evalue=1e-10,
        )
        self.assertEqual([(hit.qseqid, hit.sseqid) for hit in filtered], [("P1", "Q1")])

    def test_rbh_is_selected_at_gene_level(self) -> None:
        mapping1 = {
            "P1": homolog.ProteinMapping("A", "a"),
            "P2": homolog.ProteinMapping("A", "a-isoform"),
            "P3": homolog.ProteinMapping("B", "b"),
        }
        mapping2 = {
            "Q1": homolog.ProteinMapping("X", "x"),
            "Q2": homolog.ProteinMapping("Y", "y"),
        }
        forward = [
            homolog.Hit("P1", "Q1", 80, 100, 100, 100, 100, 100, 1e-30, 100),
            homolog.Hit("P2", "Q1", 90, 100, 100, 100, 100, 100, 1e-30, 110),
            homolog.Hit("P1", "Q2", 95, 100, 100, 100, 100, 100, 1e-30, 120),
            homolog.Hit("P3", "Q2", 90, 100, 100, 100, 100, 100, 1e-30, 90),
        ]
        reverse = [
            homolog.Hit("Q1", "P2", 90, 100, 100, 100, 100, 100, 1e-30, 110),
            homolog.Hit("Q2", "P1", 95, 100, 100, 100, 100, 100, 1e-30, 120),
        ]
        pairs, forward_count, reverse_count = homolog.select_pairs(
            forward,
            reverse,
            mapping1,
            mapping2,
            mode="rbh",
            min_identity=40,
            min_query_cover=70,
            min_subject_cover=70,
            evalue=1e-10,
        )
        self.assertEqual((forward_count, reverse_count), (4, 2))
        self.assertEqual([(pair.gene1, pair.gene2) for pair in pairs], [("A", "Y")])
        self.assertEqual(pairs[0].hit.sseqid, "Q2")

    def test_all_mode_keeps_many_to_many_gene_pairs(self) -> None:
        mapping1 = {
            "P1": homolog.ProteinMapping("A", "a"),
            "P2": homolog.ProteinMapping("A", "a"),
        }
        mapping2 = {
            "Q1": homolog.ProteinMapping("X", "x"),
            "Q2": homolog.ProteinMapping("Y", "y"),
        }
        hits = [
            homolog.Hit("P1", "Q1", 80, 100, 100, 100, 100, 100, 1e-20, 100),
            homolog.Hit("P2", "Q2", 80, 100, 100, 100, 100, 100, 1e-20, 90),
        ]
        pairs, _, _ = homolog.select_pairs(
            hits,
            [],
            mapping1,
            mapping2,
            mode="all",
            min_identity=40,
            min_query_cover=70,
            min_subject_cover=70,
            evalue=1e-10,
        )
        self.assertEqual([(pair.gene1, pair.gene2) for pair in pairs], [("A", "X"), ("A", "Y")])

    def test_run_analysis_uses_two_searches_and_writes_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta1 = root / "one.faa"
            fasta2 = root / "two.faa"
            fasta1.write_text(">P1 gene-one [GeneID=123]\nMKTK\n", encoding="utf-8")
            fasta2.write_text(
                ">FBpp0000001 type=polypeptide; ID=FBpp0000001; parent=FBgn0000001; name=gene-two;\nMKTK\n",
                encoding="utf-8",
            )
            output = root / "pairs.tsv"
            calls: list[list[str]] = []

            def fake_run(command: list[str]) -> homolog.subprocess.CompletedProcess[str]:
                calls.append(command)
                if command[1] == "version":
                    return homolog.subprocess.CompletedProcess(
                        command, 0, "diamond version 2.2.5", ""
                    )
                if command[1] == "makedb":
                    return homolog.subprocess.CompletedProcess(command, 0, "", "")
                out = Path(command[command.index("--out") + 1])
                query = Path(command[command.index("--query") + 1]).name
                if query == fasta1.name:
                    line = hit_line("P1", "FBpp0000001", length=4, qlen=4, slen=4, bitscore=20)
                else:
                    line = hit_line("FBpp0000001", "P1", length=4, qlen=4, slen=4, bitscore=20)
                out.write_text(line + "\n", encoding="utf-8")
                return homolog.subprocess.CompletedProcess(command, 0, "", "")

            with patch.object(homolog, "_run_command", side_effect=fake_run):
                summary = homolog.run_analysis(
                    fasta1,
                    fasta2,
                    output=output,
                    mode="rbh",
                    min_identity=40,
                    min_query_cover=70,
                    min_subject_cover=70,
                    evalue=1e-10,
                    threads=1,
                    diamond="fake-diamond",
                )
            self.assertEqual(summary.pair_count, 1)
            self.assertEqual(len(calls), 5)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["gene1"], "123")
            self.assertEqual(rows[0]["gene2"], "FBgn0000001")
            self.assertEqual(rows[0]["homology"], "80")
            metadata = json.loads(summary.metadata.read_text(encoding="utf-8"))
            self.assertEqual(metadata["counts"]["selected_pairs"], 1)
            self.assertEqual(metadata["diamond"]["version"], "diamond version 2.2.5")

    def test_unmapped_proteins_fail_unless_explicitly_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "one.faa"
            path.write_text(">P1\nMKT\n", encoding="utf-8")
            with self.assertRaisesRegex(homolog.HomologyError, "Could not map"):
                homolog._resolve_mapping(
                    path,
                    homolog._read_fasta(path),
                    manifest=None,
                    explicit_map=None,
                    allow_unmapped=False,
                )

    def test_empty_diamond_output_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.tsv"
            path.write_text("", encoding="utf-8")
            self.assertEqual(homolog._parse_hits(path), [])

    def test_duplicate_fasta_accessions_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "duplicate.faa"
            path.write_text(">P1\nMKT\n>P1\nMKTK\n", encoding="utf-8")
            with self.assertRaisesRegex(homolog.HomologyError, "Duplicate protein"):
                homolog._read_fasta(path)

    def test_diamond_failure_is_reported(self) -> None:
        result = homolog.subprocess.CompletedProcess(
            ["diamond", "blastp"], 2, "", "invalid database"
        )
        with self.assertRaisesRegex(homolog.HomologyError, "invalid database"):
            homolog._check_command(result, list(result.args))


if __name__ == "__main__":
    unittest.main()
