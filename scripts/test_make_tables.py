import csv
import tempfile
import unittest
from pathlib import Path

from scripts import make_tables


ROOT = Path(__file__).resolve().parents[1]
ROUND_SOURCE = ROOT / "round-numbers-overview.txt"
MANIFEST = ROOT / "benchmark-manifest.csv"
# Independent of make_tables.EXPECTED_REPETITIONS (100 in production) so the
# fixture stays small and tests stay fast; the validation logic is generic
# over the repetition count either way.
TEST_REPETITIONS = 5


class MakeTablesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = make_tables.parse_round_source(ROUND_SOURCE)
        cls.expected = make_tables.expected_cases(cls.source)

    def write_csv(self, directory, mutate=None):
        rows = []
        for case_index, (key, expected) in enumerate(self.expected.items()):
            for rep in range(TEST_REPETITIONS):
                rows.append(
                    {
                        "construction": key.construction,
                        "variant": key.variant,
                        "field": key.field,
                        "t": "-" if key.t is None else str(key.t),
                        "round_set": key.round_set,
                        "rounds": expected.rounds,
                        "iters": "10",
                        "rep": str(rep),
                        "total_ns": str(10 * (2000 + case_index * 10 + rep)),
                    }
                )
        if mutate is not None:
            mutate(rows)
        path = Path(directory) / "results.csv"
        with path.open("w", encoding="utf-8", newline="") as output:
            output.write("# cpu: test fixture\n")
            writer = csv.DictWriter(output, fieldnames=make_tables.CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_source_derived_cases_match_checked_in_manifest(self):
        with MANIFEST.open(encoding="utf-8", newline="") as manifest_file:
            reader = csv.DictReader(
                line for line in manifest_file if line.strip() and not line.startswith("#")
            )
            manifest = {
                (
                    row["construction"],
                    row["variant"],
                    row["field"],
                    row["t"],
                    row["round_set"],
                    row["rounds"],
                )
                for row in reader
            }
        derived = {
            (
                key.construction,
                key.variant,
                key.field,
                "-" if key.t is None else str(key.t),
                key.round_set,
                expected.rounds,
            )
            for key, expected in self.expected.items()
        }
        self.assertEqual(manifest, derived)

    def test_valid_csv_generates_four_tables_and_pair_suffixes(self):
        with tempfile.TemporaryDirectory() as directory:
            input_path = self.write_csv(directory)
            loaded = make_tables.load_csv(
                input_path, self.expected, expected_repetitions=TEST_REPETITIONS
            )
            self.assertEqual([], loaded.issues)
            latex = make_tables.generate_document(
                input_path, ROUND_SOURCE, self.source, loaded
            )
        self.assertEqual(4, latex.count(r"\begin{table}[htb]"))
        self.assertNotIn(r"\textit{Type", latex)
        self.assertIn(r"\begin{tabular}{l|S", latex)
        self.assertIn(r"@{}l|S", latex)
        self.assertIn(r"\multicolumn{4}{c}{BLS12-381}", latex)
        self.assertIn(r"\multicolumn{2}{c}{\textemdash}", latex)
        self.assertNotIn("|c}", latex)
        self.assertNotIn("||", latex)
        self.assertIn(r"\xhasheight", latex)
        self.assertIn(r"\xhashtwentyfour", latex)
        self.assertIn(r"\,(\num{", latex)
        self.assertIn(
            "the original round count data is given in parentheses.", latex
        )
        self.assertNotIn("runtime at the designers' original round count", latex)
        self.assertIn("Benchmark provenance: cpu: test fixture", latex)
        self.assertRegex(latex, r"\$t=3\$ \(\\si\{\\(nano|micro)\\second\}\)")

    def test_rows_are_sorted_slowest_first(self):
        def mutate(rows):
            for row in rows:
                if (
                    row["construction"] == "poseidon"
                    and row["field"] == "bls12_381"
                    and row["t"] == "3"
                ):
                    row["total_ns"] = str(10 * (9_000_000 + int(row["rep"])))

        with tempfile.TemporaryDirectory() as directory:
            input_path = self.write_csv(directory, mutate)
            loaded = make_tables.load_csv(
                input_path, self.expected, expected_repetitions=TEST_REPETITIONS
            )
            latex = make_tables.generate_document(
                input_path, ROUND_SOURCE, self.source, loaded
            )
        first_table = latex.split(r"\end{table}", 1)[0]
        self.assertLess(
            first_table.index(r"\poseidon &"), first_table.index(r"\gmimchash &")
        )

    def test_excluded_monolith_cells_are_empty_without_warnings(self):
        with tempfile.TemporaryDirectory() as directory:
            loaded = make_tables.load_csv(
                self.write_csv(directory), self.expected, expected_repetitions=TEST_REPETITIONS
            )
        self.assertFalse(any("monolith/-/koalabear" in issue for issue in loaded.issues))
        cell = self.source[("monolith", "koalabear", 16)]
        rendered = make_tables.latex_cell(cell, "", loaded.measurements, 1000.0, 2)
        self.assertEqual(r"\multicolumn{2}{c}{}", rendered)

    def test_missing_duplicate_wrong_round_and_extra_cases_are_reported(self):
        def mutate(rows):
            rows.pop(0)
            rows[0]["rounds"] = "999"
            rows.append(dict(rows[1]))
            extra = dict(rows[2])
            extra["construction"] = "not_in_manifest"
            rows.append(extra)

        with tempfile.TemporaryDirectory() as directory:
            loaded = make_tables.load_csv(
                self.write_csv(directory, mutate),
                self.expected,
                expected_repetitions=TEST_REPETITIONS,
            )
        issues = "\n".join(loaded.issues)
        self.assertIn("wrong rounds", issues)
        self.assertIn("duplicate rep", issues)
        self.assertIn("extra CSV case", issues)
        self.assertIn("has reps", issues)

    def test_unit_boundary_is_strict_and_precision_avoids_collisions(self):
        self.assertEqual("ns", make_tables.numeric_format([1000.0])[0])
        self.assertEqual("us", make_tables.numeric_format([1000.0001])[0])
        precision = make_tables.choose_precision([1000.001, 1000.004], 1000.0)
        self.assertGreater(precision, 2)


if __name__ == "__main__":
    unittest.main()
