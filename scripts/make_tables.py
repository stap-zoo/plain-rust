#!/usr/bin/env python3
"""Validate publication benchmark CSV and generate the four LaTeX tables."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_REPETITIONS = 100
CSV_COLUMNS = (
    "construction",
    "variant",
    "field",
    "t",
    "round_set",
    "rounds",
    "iters",
    "rep",
    "total_ns",
)
FIELDS_BLS_BN254 = (
    ("bls12_381", 3, "BLS12-381", "3"),
    ("bls12_381", 4, "BLS12-381", "4"),
    ("bn254", 3, "BN254", "3"),
    ("bn254", 4, "BN254", "4"),
)
FIELDS_GOLDILOCKS_MERSENNE31 = (
    ("goldilocks", 8, "Goldilocks", "8"),
    ("goldilocks", 12, "Goldilocks", "12"),
    ("mersenne31", 16, "Mersenne31", "16"),
    ("mersenne31", 24, "Mersenne31", "24"),
)
FIELDS_KOALABEAR_BABYBEAR = (
    ("koalabear", 16, "KoalaBear", "16"),
    ("koalabear", 24, "KoalaBear", "24"),
    ("babybear", 16, "BabyBear", "16"),
    ("babybear", 24, "BabyBear", "24"),
)
ALL_FIELDS = FIELDS_BLS_BN254 + FIELDS_GOLDILOCKS_MERSENNE31 + FIELDS_KOALABEAR_BABYBEAR

CONSTRUCTION_KEYS = {
    "GMiMCHash": "gmimc",
    "GMiMCHash2": "gmimc2",
    "Poseidon": "poseidon",
    "Poseidon2": "poseidon2",
    "Neptune": "neptune",
    "Rescue-Prime": "rescue_prime",
    "Anemoi": "anemoi",
    "Griffin": "griffin",
    "Arion": "arion",
    "Grendel": "grendel",
    "Reinforced Concrete": "reinforced_concrete",
    "Skyscraper": "skyscraper",
    "Monolith": "monolith",
    "Polocolo": "polocolo",
    "XHash8/12": "xhash",
    "Tip4'": "tip4p",
    "pSquare-Hash": "psquarehash",
}
DISPLAY_NAMES = {
    "gmimc": r"\gmimchash",
    "gmimc2": r"\gmimchashtwo",
    "poseidon": r"\poseidon",
    "poseidon2": r"\poseidontwo",
    "neptune": r"\neptune",
    "psquarehash": r"\psquarehash",
    "rescue_prime": r"\rescueprime",
    "anemoi": r"\anemoi",
    "griffin": r"\griffin",
    "arion": r"\arion",
    "grendel": r"\grendel",
    "xhash": r"\xhash",
    "reinforced_concrete": r"\reinforcedc",
    "skyscraper": r"\skyscraper",
    "monolith": r"\monolith",
    "polocolo": r"\polocolo",
    "tip4p": r"\tipfourprime",
}
VARIANT_DISPLAY_NAMES = {
    "XHash8": r"\xhasheight",
    "XHash12": r"\xhashtwelve",
    "XHash16": r"\xhashsixteen",
    "XHash24": r"\xhashtwentyfour",
}
CONSTRUCTIONS = (
    "gmimc",
    "gmimc2",
    "poseidon",
    "poseidon2",
    "neptune",
    "psquarehash",
    "rescue_prime",
    "anemoi",
    "griffin",
    "arion",
    "grendel",
    "xhash",
    "reinforced_concrete",
    "skyscraper",
    "monolith",
    "polocolo",
    "tip4p",
)
PLAIN = (
    ("sha256_compress", "SHA-256 compression", "256", "512"),
    ("keccak_f1600", r"Keccak-$f[1600]$", "1600", r"\textemdash"),
    ("blake2b_compress", "BLAKE2b compression", "512", "1024"),
    ("blake3_compress", "BLAKE3 compression", "256", "512"),
)

# Section 2's sole defined-but-unimplemented cells. They intentionally render empty.
EXCLUDED = {
    ("monolith", "koalabear", 16),
    ("monolith", "koalabear", 24),
    ("monolith", "babybear", 16),
    ("monolith", "babybear", 24),
}


@dataclass(frozen=True, order=True)
class CaseKey:
    construction: str
    variant: str
    field: str
    t: int | None
    round_set: str

    def describe(self) -> str:
        width = "-" if self.t is None else str(self.t)
        return "/".join(
            (self.construction, self.variant or "-", self.field, width, self.round_set)
        )


@dataclass(frozen=True)
class SourceCell:
    construction: str
    field: str
    t: int
    current: str | None
    original: str | None = None

    @property
    def undefined(self) -> bool:
        return self.current is None


@dataclass(frozen=True)
class Measurement:
    median_ns: float
    min_ns: float
    max_ns: float
    spread: float


@dataclass(frozen=True)
class Expected:
    rounds: str


@dataclass
class LoadedCsv:
    measurements: dict[CaseKey, Measurement]
    rounds: dict[CaseKey, str]
    preamble: list[str]
    issues: list[str]


def normalize_construction(raw: str) -> str:
    name = raw.replace("‡", "").strip()
    name = name.replace("’", "'").replace("′", "'")
    if name not in CONSTRUCTION_KEYS:
        raise ValueError(f"unknown construction in round source: {raw!r}")
    return CONSTRUCTION_KEYS[name]


def normalize_rounds(raw: str) -> str:
    return re.sub(r"\s*\+\s*", "+", raw.strip())


def parse_source_value(construction: str, raw: str) -> tuple[str | None, str | None]:
    value = raw.strip()
    if value in {"—", "--", "---", r"\textemdash"}:
        return None, None
    if construction == "gmimc2":
        match = re.fullmatch(r"(\d+)\s*\([^)]*\)", value)
        if not match:
            raise ValueError(f"invalid GMiMCHash2 source cell: {raw!r}")
        return match.group(1), None
    match = re.fullmatch(r"(.+?)\s*\(([^()]*)\)", value)
    if match:
        return normalize_rounds(match.group(1)), normalize_rounds(match.group(2))
    value = normalize_rounds(value)
    if not re.fullmatch(r"\d+(?:\+\d+)?", value):
        raise ValueError(f"invalid round source cell: {raw!r}")
    return value, None


def parse_header_cell(raw: str) -> tuple[str, int]:
    match = re.fullmatch(r"(.+?)\s+\$t=(\d+)\$", raw.strip())
    if not match:
        raise ValueError(f"invalid round-source column header: {raw!r}")
    field_name = match.group(1).strip()
    field = {
        "BLS12-381": "bls12_381",
        "BN254": "bn254",
        "Goldilocks": "goldilocks",
        "Mersenne31": "mersenne31",
        "KoalaBear": "koalabear",
        "BabyBear": "babybear",
    }.get(field_name)
    if field is None:
        raise ValueError(f"unknown field in round source: {field_name!r}")
    return field, int(match.group(2))


def markdown_cells(line: str) -> list[str]:
    return [part.strip() for part in line.strip().strip("|").split("|")]


def parse_round_source(path: Path) -> dict[tuple[str, str, int], SourceCell]:
    lines = path.read_text(encoding="utf-8").splitlines()
    cells: dict[tuple[str, str, int], SourceCell] = {}
    index = 0
    tables = 0
    while index < len(lines):
        if not lines[index].startswith("| Construction |"):
            index += 1
            continue
        headers = markdown_cells(lines[index])
        columns = [parse_header_cell(value) for value in headers[1:]]
        index += 2  # Skip the Markdown separator.
        tables += 1
        while index < len(lines) and lines[index].lstrip().startswith("|"):
            row = markdown_cells(lines[index])
            if len(row) != len(headers):
                raise ValueError(
                    f"round-source row has {len(row)} columns; expected {len(headers)}: {lines[index]!r}"
                )
            construction = normalize_construction(row[0])
            for (field, width), raw in zip(columns, row[1:]):
                current, original = parse_source_value(construction, raw)
                key = (construction, field, width)
                if key in cells:
                    raise ValueError(f"duplicate round-source cell: {key}")
                cells[key] = SourceCell(construction, field, width, current, original)
            index += 1
    if tables != 2:
        raise ValueError(f"expected two round-number tables, found {tables}")
    expected_count = len(CONSTRUCTION_KEYS) * len(ALL_FIELDS)
    if len(cells) != expected_count:
        raise ValueError(f"expected {expected_count} round-source cells, found {len(cells)}")
    return cells


def cell_variants(cell: SourceCell) -> tuple[str, ...]:
    if cell.construction == "gmimc2":
        if cell.field in {"bls12_381", "bn254"}:
            return ("alpha=8",)
        if cell.field == "goldilocks":
            return ("alpha=4",)
        return ("alpha=2",)
    if cell.construction == "grendel":
        return ("native",)
    if cell.construction == "xhash":
        if (cell.field, cell.t) == ("goldilocks", 12):
            return ("XHash8", "XHash12")
        if (cell.field, cell.t) == ("mersenne31", 24):
            return ("XHash16", "XHash24")
    return ("",)


def expected_cases(
    source: dict[tuple[str, str, int], SourceCell]
) -> dict[CaseKey, Expected]:
    expected: dict[CaseKey, Expected] = {}
    for cell_key, cell in source.items():
        if cell.undefined or cell_key in EXCLUDED:
            continue
        for variant in cell_variants(cell):
            if cell.original is None:
                expected[CaseKey(cell.construction, variant, cell.field, cell.t, "single")] = Expected(
                    cell.current or ""
                )
            else:
                expected[CaseKey(cell.construction, variant, cell.field, cell.t, "current")] = Expected(
                    cell.current or ""
                )
                expected[CaseKey(cell.construction, variant, cell.field, cell.t, "original")] = Expected(
                    cell.original
                )
    for key, _display, _state, _block in PLAIN:
        expected[CaseKey("plain", key, "-", None, "single")] = Expected("-")
    return expected


def parse_optional_t(raw: str) -> int | None:
    if raw == "-":
        return None
    value = int(raw)
    if value <= 0:
        raise ValueError("t must be positive")
    return value


def load_csv(
    path: Path,
    expected: dict[CaseKey, Expected],
    expected_repetitions: int = EXPECTED_REPETITIONS,
) -> LoadedCsv:
    preamble: list[str] = []
    data_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# "):
            preamble.append(line[2:])
        elif line.startswith("#") or not line.strip():
            continue
        else:
            data_lines.append(line)
    issues: list[str] = []
    if not data_lines:
        return LoadedCsv({}, {}, preamble, ["benchmark CSV contains no data"])
    reader = csv.DictReader(data_lines)
    if reader.fieldnames != list(CSV_COLUMNS):
        return LoadedCsv(
            {},
            {},
            preamble,
            [f"CSV header must be exactly {','.join(CSV_COLUMNS)}"],
        )

    samples: dict[CaseKey, dict[int, float]] = {}
    rounds_by_key: dict[CaseKey, str] = {}
    for line_number, row in enumerate(reader, start=2):
        try:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("row has the wrong number of columns")
            key = CaseKey(
                row["construction"],
                row["variant"],
                row["field"],
                parse_optional_t(row["t"]),
                row["round_set"],
            )
            iters = int(row["iters"])
            rep = int(row["rep"])
            total_ns = int(row["total_ns"])
            if iters <= 0:
                raise ValueError("iters must be positive")
            if total_ns < 0:
                raise ValueError("total_ns must be non-negative")
            timing = total_ns / iters
            if not math.isfinite(timing):
                raise ValueError("timing is not finite")
        except (KeyError, OverflowError, TypeError, ValueError) as error:
            issues.append(f"CSV row {line_number}: {error}")
            continue
        rounds = row["rounds"]
        if key not in expected:
            issues.append(f"extra CSV case: {key.describe()}")
            continue
        wanted = expected[key].rounds
        if rounds != wanted:
            issues.append(
                f"wrong rounds for {key.describe()}: got {rounds!r}, expected {wanted!r}"
            )
            continue
        old_rounds = rounds_by_key.setdefault(key, rounds)
        if old_rounds != rounds:
            issues.append(f"inconsistent rounds for {key.describe()}")
            continue
        repetitions = samples.setdefault(key, {})
        if rep in repetitions:
            issues.append(f"duplicate rep {rep} for {key.describe()}")
            continue
        repetitions[rep] = timing

    measurements: dict[CaseKey, Measurement] = {}
    for key, expected_case in expected.items():
        repetitions = samples.get(key)
        if repetitions is None:
            issues.append(f"missing CSV case: {key.describe()}")
            continue
        actual_reps = set(repetitions)
        required_reps = set(range(expected_repetitions))
        if actual_reps != required_reps:
            issues.append(
                f"{key.describe()} has reps {sorted(actual_reps)}; "
                f"expected {expected_repetitions} reps [0, ..., {expected_repetitions - 1}]"
            )
            continue
        values = list(repetitions.values())
        median = statistics.median(values)
        minimum = min(values)
        maximum = max(values)
        spread = 0.0 if median == 0.0 else (maximum - minimum) / median
        measurements[key] = Measurement(median, minimum, maximum, spread)
        rounds_by_key.setdefault(key, expected_case.rounds)
    return LoadedCsv(measurements, rounds_by_key, preamble, issues)


def timings_for_cell(
    construction: str,
    variant: str,
    field: str,
    width: int,
    measurements: dict[CaseKey, Measurement],
) -> tuple[Measurement | None, Measurement | None]:
    single = measurements.get(CaseKey(construction, variant, field, width, "single"))
    if single is not None:
        return single, None
    current = measurements.get(CaseKey(construction, variant, field, width, "current"))
    original = measurements.get(CaseKey(construction, variant, field, width, "original"))
    return current, original


def table_variants(
    construction: str,
    fields: Sequence[tuple[str, int, str, str]],
    source: dict[tuple[str, str, int], SourceCell],
) -> tuple[str, ...]:
    if construction != "xhash":
        return ("",)
    variants: list[str] = []
    for field, width, _display, _t in fields:
        cell = source[(construction, field, width)]
        if cell.undefined:
            continue
        for variant in cell_variants(cell):
            if variant not in variants:
                variants.append(variant)
    return tuple(variants) if variants else ("",)


def row_is_relevant(
    construction: str,
    fields: Sequence[tuple[str, int, str, str]],
    source: dict[tuple[str, str, int], SourceCell],
) -> bool:
    return any(not source[(construction, field, width)].undefined for field, width, _, _ in fields)


def choose_precision(values: Iterable[float], scale: float) -> int:
    distinct = sorted(set(values))
    for precision in range(2, 13):
        rendered = [f"{value / scale:.{precision}f}" for value in distinct]
        if len(rendered) == len(set(rendered)):
            return precision
    return 12


def numeric_format(values: Sequence[float]) -> tuple[str, float, int, str]:
    microseconds = bool(values) and all(value > 1000.0 for value in values)
    scale = 1000.0 if microseconds else 1.0
    precision = choose_precision(values, scale)
    unit = r"\si{\micro\second}" if microseconds else r"\si{\nano\second}"
    return ("us" if microseconds else "ns"), scale, precision, unit


def format_number(value: float, scale: float, precision: int) -> str:
    return f"{value / scale:.{precision}f}"


def latex_cell(
    cell: SourceCell,
    variant: str,
    measurements: dict[CaseKey, Measurement],
    scale: float,
    precision: int,
) -> str:
    if cell.undefined:
        return r"\multicolumn{2}{c}{\textemdash}"
    if (cell.construction, cell.field, cell.t) in EXCLUDED:
        return r"\multicolumn{2}{c}{}"
    applicable_variants = cell_variants(cell)
    if variant and variant not in applicable_variants:
        # This row's variant belongs to a different field's special case
        # (e.g. XHash8/XHash12 are Goldilocks-only, XHash16/XHash24 are
        # Mersenne31-only); this cell simply doesn't apply to that variant.
        return r"\multicolumn{2}{c}{\textemdash}"
    actual_variant = variant or applicable_variants[0]
    primary, original = timings_for_cell(
        cell.construction, actual_variant, cell.field, cell.t, measurements
    )
    if primary is None:
        return r"\multicolumn{2}{c}{\textbf{?}}"
    primary_text = format_number(primary.median_ns, scale, precision)
    suffix = ""
    if cell.original is not None:
        if original is None:
            return r"\multicolumn{2}{c}{\textbf{?}}"
        suffix = rf"\,(\num{{{format_number(original.median_ns, scale, precision)}}})"
    return f"{primary_text} & {suffix}"


def si_setup() -> str:
    return r"\sisetup{group-separator={\,}, group-minimum-digits=4}"


def result_column_spec(precision: int, count: int) -> str:
    # Seven integer digits covers the expected native timings while preserving alignment.
    one = f"S[table-format=7.{precision}]@{{}}l"
    return "l" + one * count


def grouped_result_column_spec(
    fields: Sequence[tuple[str, int, str, str]],
    column_formats: Sequence[tuple[str, float, int, str]],
) -> str:
    parts = ["l|"]
    previous_prime = fields[0][2]
    for index, ((_field, _width, prime, _t), (_unit_kind, _scale, precision, _unit)) in enumerate(
        zip(fields, column_formats)
    ):
        if index and prime != previous_prime:
            parts.append("|")
        parts.append(f"S[table-format=7.{precision}]@{{}}l")
        previous_prime = prime
    return "".join(parts)


def sorted_table_rows(
    fields: Sequence[tuple[str, int, str, str]],
    source: dict[tuple[str, str, int], SourceCell],
    measurements: dict[CaseKey, Measurement],
) -> list[tuple[str, str]]:
    rows = [
        (construction, variant)
        for construction in CONSTRUCTIONS
        if row_is_relevant(construction, fields, source)
        for variant in table_variants(construction, fields, source)
    ]
    first_field, first_width, _display, _t = fields[0]

    def first_column_timing(row: tuple[str, str]) -> tuple[int, float]:
        construction, variant = row
        cell = source[(construction, first_field, first_width)]
        actual_variant = variant or cell_variants(cell)[0]
        primary, _original = timings_for_cell(
            construction, actual_variant, first_field, first_width, measurements
        )
        # Rows with no measurement in the first data column stay last,
        # regardless of sort direction; defined rows sort slowest-first.
        return (1, 0.0) if primary is None else (0, -primary.median_ns)

    # Python's stable sort retains the source construction order for rows that
    # have no defined measurement in the first data column.
    rows.sort(key=first_column_timing)
    return rows


def column_values(
    field: str,
    width: int,
    rows: Sequence[tuple[str, str]],
    source: dict[tuple[str, str, int], SourceCell],
    measurements: dict[CaseKey, Measurement],
) -> list[float]:
    values: list[float] = []
    for construction, variant in rows:
        cell = source[(construction, field, width)]
        if cell.undefined or (construction, field, width) in EXCLUDED:
            continue
        applicable_variants = cell_variants(cell)
        if variant and variant not in applicable_variants:
            continue
        actual_variant = variant or applicable_variants[0]
        primary, original = timings_for_cell(construction, actual_variant, field, width, measurements)
        if primary is not None:
            values.append(primary.median_ns)
        if original is not None:
            values.append(original.median_ns)
    return values


def grouped_table(
    fields: Sequence[tuple[str, int, str, str]],
    source: dict[tuple[str, str, int], SourceCell],
    measurements: dict[CaseKey, Measurement],
    *,
    caption_prefix: str,
    label: str,
) -> str:
    rows = sorted_table_rows(fields, source, measurements)
    column_formats = [
        numeric_format(column_values(field, width, rows, source, measurements))
        for field, width, _display, _t in fields
    ]
    lines = [
        r"\begin{table}[htb]",
        r"\centering",
        r"\begingroup",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2pt}",
        si_setup(),
        rf"\caption{{{caption_prefix} Medians of {EXPECTED_REPETITIONS} repetitions; units are given per column. For cells whose round count was updated, the original round count data is given in parentheses.}}",
        rf"\label{{{label}}}",
        rf"\begin{{tabular}}{{{grouped_result_column_spec(fields, column_formats)}}}",
        r"\toprule",
    ]
    field_names: list[str] = []
    for _field, _width, display, _t in fields:
        if display not in field_names:
            field_names.append(display)
    first_header = [""]
    column = 2
    cmidrules: list[str] = []
    for display in field_names:
        count = sum(1 for item in fields if item[2] == display)
        span = count * 2
        first_header.append(rf"\multicolumn{{{span}}}{{c}}{{{display}}}")
        cmidrules.append(rf"\cmidrule(lr){{{column}-{column + span - 1}}}")
        column += span
    lines.append(" & ".join(first_header) + r" \\")
    lines.append("".join(cmidrules))
    second_header = ["Construction"]
    for (_field, _width, _prime, width), (_unit_kind, _scale, _precision, unit) in zip(
        fields, column_formats
    ):
        second_header.append(rf"\multicolumn{{2}}{{c}}{{$t={width}$ ({unit})}}")
    lines.append(" & ".join(second_header) + r" \\")
    lines.append(r"\midrule")

    for construction, variant in rows:
        label_text = VARIANT_DISPLAY_NAMES.get(variant, variant) if variant else DISPLAY_NAMES[construction]
        cells = [
            latex_cell(
                source[(construction, field, width)],
                variant,
                measurements,
                scale,
                precision,
            )
            for (field, width, _prime, _t), (_unit_kind, scale, precision, _unit) in zip(
                fields, column_formats
            )
        ]
        lines.append(" & ".join([label_text] + cells) + r" \\")
    lines.extend((r"\bottomrule", r"\end{tabular}", r"\endgroup", r"\end{table}"))
    return "\n".join(lines)


def plain_table(measurements: dict[CaseKey, Measurement]) -> str:
    values = [
        measurements[key].median_ns
        for construction, _display, _state, _block in PLAIN
        if (key := CaseKey("plain", construction, "-", None, "single")) in measurements
    ]
    _unit_kind, scale, precision, unit = numeric_format(values)
    lines = [
        r"\begin{table}[htb]",
        r"\centering",
        r"\begingroup",
        si_setup(),
        rf"\caption{{Plain-hash baseline timings in {unit}, reported as medians of {EXPECTED_REPETITIONS} repetitions. State and block sizes are in bits; Keccak-$f[1600]$ is a permutation and has no fixed block size.}}",
        r"\label{tab:plain-hash-benchmarks}",
        rf"\begin{{tabular}}{{lcc{result_column_spec(precision, 1)[1:]}}}",
        r"\toprule",
        r"Construction & State & Block & \multicolumn{2}{c}{Runtime} \\",
        r"\midrule",
    ]
    for construction, display, state, block in PLAIN:
        key = CaseKey("plain", construction, "-", None, "single")
        measurement = measurements.get(key)
        timing = (
            r"\multicolumn{2}{c}{\textbf{?}}"
            if measurement is None
            else f"{format_number(measurement.median_ns, scale, precision)} & "
        )
        lines.append(f"{display} & {state} & {block} & {timing}" + r" \\")
    lines.extend((r"\bottomrule", r"\end{tabular}", r"\endgroup", r"\end{table}"))
    return "\n".join(lines)


def provenance_comments(
    input_path: Path,
    source_path: Path,
    loaded: LoadedCsv,
) -> str:
    lines = [
        "% Generated by scripts/make_tables.py; do not edit by hand.",
        f"% Benchmark input: {input_path}",
        f"% Round source: {source_path}",
        "% Requires: \\usepackage{booktabs} and \\usepackage{siunitx}",
    ]
    lines.extend(f"% Benchmark provenance: {item}" for item in loaded.preamble)
    for key in sorted(loaded.measurements):
        item = loaded.measurements[key]
        lines.append(
            f"% stats {key.describe()}: median_ns={item.median_ns:.9g}; "
            f"min_ns={item.min_ns:.9g}; max_ns={item.max_ns:.9g}; "
            f"spread={item.spread * 100:.3f}%"
        )
    return "\n".join(lines)


def generate_document(
    input_path: Path,
    source_path: Path,
    source: dict[tuple[str, str, int], SourceCell],
    loaded: LoadedCsv,
) -> str:
    sections = [
        provenance_comments(input_path, source_path, loaded),
        grouped_table(
            FIELDS_BLS_BN254,
            source,
            loaded.measurements,
            caption_prefix="Native permutation timings over BLS12-381 and BN254.",
            label="tab:benchmark-bls-bn254",
        ),
        grouped_table(
            FIELDS_GOLDILOCKS_MERSENNE31,
            source,
            loaded.measurements,
            caption_prefix="Native permutation timings over Goldilocks and Mersenne31.",
            label="tab:benchmark-goldilocks-mersenne31",
        ),
        grouped_table(
            FIELDS_KOALABEAR_BABYBEAR,
            source,
            loaded.measurements,
            caption_prefix="Native permutation timings over KoalaBear and BabyBear.",
            label="tab:benchmark-koalabear-babybear",
        ),
        plain_table(loaded.measurements),
    ]
    return "\n\n".join(sections) + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="publication benchmark CSV")
    parser.add_argument("--round-source", required=True, type=Path, help="round-number Markdown source")
    parser.add_argument("--output", required=True, type=Path, help="generated LaTeX file")
    parser.add_argument("--strict", action="store_true", help="fail instead of rendering validation gaps")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        source = parse_round_source(args.round_source)
        expected = expected_cases(source)
        loaded = load_csv(args.input, expected)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    for issue in loaded.issues:
        print(f"warning: {issue}", file=sys.stderr)
    if args.strict and loaded.issues:
        print(f"error: strict validation failed with {len(loaded.issues)} issue(s)", file=sys.stderr)
        return 2
    document = generate_document(args.input, args.round_source, source, loaded)
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(document, encoding="utf-8")
    except OSError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
