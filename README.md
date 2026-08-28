# Plain Rust hash benchmarks

This repository contains plain Rust implementations and a native benchmark harness for
ZK-friendly permutations over BLS12-381, BN254, Goldilocks, Mersenne31, KoalaBear, and BabyBear.
It also benchmarks SHA-256, Keccak-f[1600], BLAKE2b, and BLAKE3 compression/permutation baselines.

## Requirements

- A recent stable Rust toolchain
- An x86_64 CPU with BMI2/ADX (Intel Haswell/2013+, AMD Zen/2017+) — `.cargo/config.toml` unconditionally compiles in ark-ff's assembly field-arithmetic backend for these extensions, and any `cargo build`/`test`/`run` here will crash with SIGILL on older hardware
- Python 3 for validation and LaTeX table generation
- Optional: `latexmk`, `booktabs`, and `siunitx` for the standalone PDF

## Validate

```bash
cargo test --lib
python3 -m unittest -v scripts/test_make_tables.py
cargo run --release -- --list
```

## Run the benchmarks

For publication measurements, use a pinned, otherwise idle CPU core and native CPU instructions:

```bash
RUSTFLAGS='-C target-cpu=native' cargo build --release
taskset -c 2 ./target/release/sok-zk-friendly-hash-functions --csv > results.csv
```

The harness warms and calibrates every case, then records 100 repetitions. Its CSV preamble
captures the CPU, compiler, Git state, release settings, CPU affinity, and available power-state
information. Replace CPU core `2` with the isolated core selected for the benchmark machine.

For a quick structural run that is not suitable for publication:

```bash
./target/release/sok-zk-friendly-hash-functions \
    --csv --preview-iters 1 > results-preview.csv
```

## Generate tables

```bash
python3 scripts/make_tables.py \
    --input results.csv \
    --round-source round-numbers-overview.txt \
    --output tables/benchmark_tables.tex \
    --strict
```

The generated file contains the three prime-field benchmark tables and the plain-hash baseline
table. To compile the standalone preview:

```bash
(cd tables && latexmk -pdf -interaction=nonstopmode -halt-on-error \
    benchmark_tables-document.tex)
```

The benchmark registry is checked against `benchmark-manifest.csv`. The table generator checks
that all expected cases have exactly 100 repetitions and that current/original round pairs match
`round-numbers-overview.txt`.
