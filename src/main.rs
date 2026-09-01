//! Publication benchmark harness for the SoK round-number tables.

use sok_zk_friendly_hash_functions::anemoi::{anemoi::Anemoi, instances::*};
use sok_zk_friendly_hash_functions::arion::{arion::Arion, instances::*};
use sok_zk_friendly_hash_functions::fields::{FieldElement, PrimeFieldMontgomery};
use sok_zk_friendly_hash_functions::gmimc::{gmimc::Gmimc, instances::*};
use sok_zk_friendly_hash_functions::gmimc2::{gmimc2::Gmimc2, instances::*};
use sok_zk_friendly_hash_functions::grendel::{grendel::Grendel, instances::*};
use sok_zk_friendly_hash_functions::griffin::{griffin::Griffin, instances::*};
use sok_zk_friendly_hash_functions::monolith::{
    instances::*,
    monolith::{Monolith31, Monolith64},
};
use sok_zk_friendly_hash_functions::neptune::{instances::*, neptune::Neptune};
use sok_zk_friendly_hash_functions::plain_hashes;
use sok_zk_friendly_hash_functions::polocolo::{instances::*, polocolo::Polocolo};
use sok_zk_friendly_hash_functions::poseidon::{instances::*, poseidon::Poseidon};
use sok_zk_friendly_hash_functions::poseidon2::{instances::*, poseidon2::Poseidon2};
use sok_zk_friendly_hash_functions::psquarehash::*;
use sok_zk_friendly_hash_functions::reinforced_concrete::{
    instances::*, reinforced_concrete::ReinforcedConcrete,
};
use sok_zk_friendly_hash_functions::rescueprime::{instances::*, rescue_prime::RescuePrime};
use sok_zk_friendly_hash_functions::skyscraper::{instances::*, skyscraper::Skyscraper};
use sok_zk_friendly_hash_functions::tip4::{instances::TIP4P_GOLDILOCKS_PARAMS, tip4::Tip4};
use sok_zk_friendly_hash_functions::xhash::{
    xhash12_goldilocks, xhash16_m31, xhash24_m31, xhash8_goldilocks,
};
use std::collections::HashSet;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::hint::black_box;
use std::process::{self, Command};
use std::time::{Duration, Instant};

const CALIBRATION_FLOOR: usize = 64;
const CALIBRATION_TARGET: Duration = Duration::from_millis(200);
const REPETITIONS: usize = 100;
const EXPECTED_MANIFEST: &str = include_str!("../benchmark-manifest.csv");

struct Case {
    construction: &'static str,
    variant: &'static str,
    field: &'static str,
    t: Option<usize>,
    round_set: &'static str,
    rounds: &'static str,
    run: Box<dyn Fn(usize) -> u128>,
}

impl Case {
    fn metadata_csv(&self) -> String {
        format!(
            "{},{},{},{},{},{}",
            self.construction,
            self.variant,
            self.field,
            display_t(self.t),
            self.round_set,
            self.rounds
        )
    }

    fn matches(&self, filter: &str) -> bool {
        self.metadata_csv()
            .to_ascii_lowercase()
            .contains(&filter.to_ascii_lowercase())
    }
}

#[derive(Default)]
struct Options {
    help: bool,
    list: bool,
    pretty: bool,
    filters: Vec<String>,
    preview_iters: Option<usize>,
}

fn main() {
    if let Err(error) = run_main() {
        eprintln!("error: {error}");
        process::exit(2);
    }
}

fn run_main() -> Result<(), String> {
    let options = parse_options(env::args().skip(1))?;
    if options.help {
        println!("Usage: sok-zk-friendly-hash-functions [--csv|--pretty|--list] [--filter TEXT ...] [--preview-iters N]\n\nCSV is the default. Repeated filters are ANDed and match all metadata columns. --preview-iters bypasses calibration and is only for non-publication table previews.");
        return Ok(());
    }
    let cases = registry();
    validate_manifest(&cases)?;
    let selected: Vec<&Case> = cases
        .iter()
        .filter(|case| options.filters.iter().all(|filter| case.matches(filter)))
        .collect();
    if selected.is_empty() {
        return Err("no benchmark cases matched the supplied filter(s)".to_owned());
    }

    if options.list {
        println!("construction,variant,field,t,round_set,rounds");
        for case in selected {
            println!("{}", case.metadata_csv());
        }
        return Ok(());
    }

    if let Some(iters) = options.preview_iters {
        eprintln!(
            "warning: preview mode uses {iters} iteration(s) per case; results are not publication measurements"
        );
    }
    print_preamble(&options);
    if !options.pretty {
        println!("construction,variant,field,t,round_set,rounds,iters,rep,total_ns");
    }

    let mut calibrated = Vec::with_capacity(selected.len());
    for case in &selected {
        (case.run)(1); // Warm lazy state and the code path before calibration.
        calibrated.push(options.preview_iters.unwrap_or_else(|| calibrate(case)));
    }

    // Repetitions are outermost so adjacent cases see comparable conditions.
    for rep in 0..REPETITIONS {
        for (case, &iters) in selected.iter().zip(&calibrated) {
            let total_ns = (case.run)(iters);
            if options.pretty {
                println!(
                    "{}/{}/{}/{}/{}/{}: rep {rep}, {iters} iters, {total_ns} ns total",
                    case.construction,
                    if case.variant.is_empty() {
                        "-"
                    } else {
                        case.variant
                    },
                    case.field,
                    display_t(case.t),
                    case.round_set,
                    case.rounds,
                );
            } else {
                println!("{},{iters},{rep},{total_ns}", case.metadata_csv());
            }
        }
    }
    Ok(())
}

fn parse_options(args: impl Iterator<Item = String>) -> Result<Options, String> {
    let mut options = Options::default();
    let mut args = args.peekable();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--csv" => options.pretty = false,
            "--pretty" => options.pretty = true,
            "--list" => options.list = true,
            "--preview-iters" => {
                options.preview_iters =
                    Some(parse_preview_iters(&args.next().ok_or_else(|| {
                        "--preview-iters requires a value".to_owned()
                    })?)?);
            }
            "--filter" => options.filters.push(
                args.next()
                    .ok_or_else(|| "--filter requires a value".to_owned())?,
            ),
            "--help" | "-h" => {
                options.help = true;
            }
            _ if arg.starts_with("--filter=") => {
                options.filters.push(arg["--filter=".len()..].to_owned());
            }
            _ if arg.starts_with("--preview-iters=") => {
                options.preview_iters =
                    Some(parse_preview_iters(&arg["--preview-iters=".len()..])?);
            }
            _ => return Err(format!("unknown argument: {arg}")),
        }
    }
    Ok(options)
}

fn parse_preview_iters(raw: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .ok()
        .filter(|&iters| iters > 0)
        .ok_or_else(|| "--preview-iters must be a positive integer".to_owned())
}

fn calibrate(case: &Case) -> usize {
    let mut iters = CALIBRATION_FLOOR;
    loop {
        if (case.run)(iters) >= CALIBRATION_TARGET.as_nanos() {
            return iters;
        }
        iters = iters
            .checked_mul(2)
            .expect("calibration iteration overflow");
    }
}

fn validate_manifest(cases: &[Case]) -> Result<(), String> {
    validate_round_pair_order(cases)?;
    let expected: Vec<&str> = EXPECTED_MANIFEST
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .skip(1)
        .collect();
    let actual: Vec<String> = cases.iter().map(Case::metadata_csv).collect();
    let expected_set: HashSet<&str> = expected.iter().copied().collect();
    let actual_set: HashSet<&str> = actual.iter().map(String::as_str).collect();
    if expected_set.len() != expected.len() {
        return Err("benchmark-manifest.csv contains duplicate rows".to_owned());
    }
    if actual_set.len() != actual.len() {
        return Err("benchmark registry contains duplicate rows".to_owned());
    }
    if expected
        .iter()
        .copied()
        .eq(actual.iter().map(String::as_str))
    {
        return Ok(());
    }
    let mut message = String::from("benchmark registry does not match benchmark-manifest.csv");
    for row in expected_set.difference(&actual_set) {
        let _ = write!(message, "\n  missing: {row}");
    }
    for row in actual_set.difference(&expected_set) {
        let _ = write!(message, "\n  extra:   {row}");
    }
    if expected_set == actual_set && expected.len() == actual.len() {
        message.push_str("\n  rows contain the same cases but registry ordering differs");
    }
    Err(message)
}

fn validate_round_pair_order(cases: &[Case]) -> Result<(), String> {
    for (index, case) in cases.iter().enumerate() {
        match case.round_set {
            "single" => {}
            "current" => {
                let Some(next) = cases.get(index + 1) else {
                    return Err(format!(
                        "current row has no adjacent original row: {}",
                        case.metadata_csv()
                    ));
                };
                if next.round_set != "original"
                    || next.construction != case.construction
                    || next.variant != case.variant
                    || next.field != case.field
                    || next.t != case.t
                {
                    return Err(format!(
                        "current row is not followed by its original pair: {}",
                        case.metadata_csv()
                    ));
                }
            }
            "original" => {
                if index == 0 || cases[index - 1].round_set != "current" {
                    return Err(format!(
                        "original row is not preceded by its current pair: {}",
                        case.metadata_csv()
                    ));
                }
            }
            other => {
                return Err(format!(
                    "invalid round_set {other:?}: {}",
                    case.metadata_csv()
                ))
            }
        }
    }
    Ok(())
}

fn display_t(t: Option<usize>) -> String {
    t.map_or_else(|| "-".to_owned(), |value| value.to_string())
}

fn add_permutation<F, R, P>(
    cases: &mut Vec<Case>,
    metadata: (
        &'static str,
        &'static str,
        &'static str,
        usize,
        &'static str,
        &'static str,
    ),
    permutation: P,
) where
    F: FieldElement + 'static,
    R: 'static,
    P: Fn(&[F]) -> R + 'static,
{
    let (construction, variant, field, t, round_set, rounds) = metadata;
    let input = make_input::<F>(t);
    let run = move |iters: usize| {
        let start = Instant::now();
        for _ in 0..iters {
            black_box(permutation(black_box(input.as_slice())));
        }
        start.elapsed().as_nanos()
    };
    cases.push(Case {
        construction,
        variant,
        field,
        t: Some(t),
        round_set,
        rounds,
        run: Box::new(run),
    });
}

macro_rules! add_case {
    ($cases:expr, $construction:literal, $variant:literal, $field:literal, $t:literal,
     $round_set:literal, $rounds:literal, $permutation:expr) => {{
        let permutation = $permutation;
        assert_eq!(permutation.get_t(), $t);
        add_permutation(
            $cases,
            ($construction, $variant, $field, $t, $round_set, $rounds),
            move |input| permutation.permutation(input),
        );
    }};
}

fn make_input<F: FieldElement>(t: usize) -> Vec<F> {
    (1..=t).map(|value| F::from_u64(value as u64)).collect()
}

fn registry() -> Vec<Case> {
    let mut cases = Vec::new();
    add_type1_cases(&mut cases);
    add_type2_cases(&mut cases);
    add_type3_cases(&mut cases);
    add_plain_cases(&mut cases);
    cases
}

fn add_type1_cases(cases: &mut Vec<Case>) {
    add_case!(
        cases,
        "gmimc",
        "",
        "bls12_381",
        3,
        "single",
        "228",
        Gmimc::new(&GMIMC_BLS12_381_3_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "bls12_381",
        4,
        "single",
        "231",
        Gmimc::new(&GMIMC_BLS12_381_4_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "bn254",
        3,
        "single",
        "228",
        Gmimc::new(&GMIMC_BN254_3_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "bn254",
        4,
        "single",
        "231",
        Gmimc::new(&GMIMC_BN254_4_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "goldilocks",
        8,
        "single",
        "68",
        Gmimc::new(&GMIMC_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "goldilocks",
        12,
        "single",
        "93",
        Gmimc::new(&GMIMC_GOLDILOCKS_12_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "mersenne31",
        16,
        "single",
        "158",
        Gmimc::new(&GMIMC_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "mersenne31",
        24,
        "single",
        "335",
        Gmimc::new(&GMIMC_MERSENNE31_24_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "koalabear",
        16,
        "single",
        "158",
        Gmimc::new(&GMIMC_KOALABEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "koalabear",
        24,
        "single",
        "335",
        Gmimc::new(&GMIMC_KOALABEAR_24_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "babybear",
        16,
        "single",
        "158",
        Gmimc::new(&GMIMC_BABYBEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "gmimc",
        "",
        "babybear",
        24,
        "single",
        "335",
        Gmimc::new(&GMIMC_BABYBEAR_24_PARAMS)
    );

    add_case!(
        cases,
        "gmimc2",
        "alpha=8",
        "bls12_381",
        3,
        "single",
        "63",
        Gmimc2::new(&GMIMC2_BLS12_381_3_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=8",
        "bls12_381",
        4,
        "single",
        "64",
        Gmimc2::new(&GMIMC2_BLS12_381_4_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=8",
        "bn254",
        3,
        "single",
        "63",
        Gmimc2::new(&GMIMC2_BN254_3_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=8",
        "bn254",
        4,
        "single",
        "64",
        Gmimc2::new(&GMIMC2_BN254_4_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=4",
        "goldilocks",
        8,
        "single",
        "88",
        Gmimc2::new(&GMIMC2_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=4",
        "goldilocks",
        12,
        "single",
        "96",
        Gmimc2::new(&GMIMC2_GOLDILOCKS_12_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=2",
        "mersenne31",
        16,
        "single",
        "176",
        Gmimc2::new(&GMIMC2_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=2",
        "mersenne31",
        24,
        "single",
        "264",
        Gmimc2::new(&GMIMC2_MERSENNE31_24_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=2",
        "koalabear",
        16,
        "single",
        "176",
        Gmimc2::new(&GMIMC2_KOALABEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=2",
        "koalabear",
        24,
        "single",
        "264",
        Gmimc2::new(&GMIMC2_KOALABEAR_24_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=2",
        "babybear",
        16,
        "single",
        "176",
        Gmimc2::new(&GMIMC2_BABYBEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "gmimc2",
        "alpha=2",
        "babybear",
        24,
        "single",
        "264",
        Gmimc2::new(&GMIMC2_BABYBEAR_24_PARAMS)
    );

    add_poseidon_cases(cases);
    add_poseidon2_cases(cases);

    add_case!(
        cases,
        "neptune",
        "",
        "bls12_381",
        4,
        "single",
        "6+69",
        Neptune::new(&NEPTUNE_BLS12_381_4_PARAMS)
    );
    add_case!(
        cases,
        "neptune",
        "",
        "bn254",
        4,
        "single",
        "6+69",
        Neptune::new(&NEPTUNE_BN254_4_PARAMS)
    );
    add_case!(
        cases,
        "neptune",
        "",
        "goldilocks",
        8,
        "single",
        "6+38",
        Neptune::new(&NEPTUNE_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "neptune",
        "",
        "goldilocks",
        12,
        "single",
        "6+42",
        Neptune::new(&NEPTUNE_GOLDILOCKS_12_PARAMS)
    );

    add_case!(
        cases,
        "psquarehash",
        "",
        "mersenne31",
        16,
        "single",
        "52",
        PSquareHash::new(&PSQUAREHASH_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "psquarehash",
        "",
        "mersenne31",
        24,
        "single",
        "52",
        PSquareHash::new(&PSQUAREHASH_MERSENNE31_24_PARAMS)
    );
    add_case!(
        cases,
        "psquarehash",
        "",
        "koalabear",
        16,
        "single",
        "52",
        PSquareHash::new(&PSQUAREHASH_KOALABEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "psquarehash",
        "",
        "koalabear",
        24,
        "single",
        "52",
        PSquareHash::new(&PSQUAREHASH_KOALABEAR_24_PARAMS)
    );
    add_case!(
        cases,
        "psquarehash",
        "",
        "babybear",
        16,
        "single",
        "52",
        PSquareHash::new(&PSQUAREHASH_BABYBEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "psquarehash",
        "",
        "babybear",
        24,
        "single",
        "52",
        PSquareHash::new(&PSQUAREHASH_BABYBEAR_24_PARAMS)
    );
}

fn add_poseidon_cases(cases: &mut Vec<Case>) {
    add_case!(
        cases,
        "poseidon",
        "",
        "bls12_381",
        3,
        "single",
        "8+56",
        Poseidon::new(&POSEIDON_BLS12_381_3_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "bls12_381",
        4,
        "single",
        "8+56",
        Poseidon::new(&POSEIDON_BLS12_381_4_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "bn254",
        3,
        "single",
        "8+56",
        Poseidon::new(&POSEIDON_BN254_3_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "bn254",
        4,
        "single",
        "8+56",
        Poseidon::new(&POSEIDON_BN254_4_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "goldilocks",
        8,
        "single",
        "8+22",
        Poseidon::new(&POSEIDON_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "goldilocks",
        12,
        "single",
        "8+22",
        Poseidon::new(&POSEIDON_GOLDILOCKS_12_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "mersenne31",
        16,
        "single",
        "8+14",
        Poseidon::new(&POSEIDON_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "mersenne31",
        24,
        "single",
        "8+22",
        Poseidon::new(&POSEIDON_MERSENNE31_24_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "koalabear",
        16,
        "single",
        "8+20",
        Poseidon::new(&POSEIDON_KOALABEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "koalabear",
        24,
        "single",
        "8+23",
        Poseidon::new(&POSEIDON_KOALABEAR_24_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "babybear",
        16,
        "single",
        "8+13",
        Poseidon::new(&POSEIDON_BABYBEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "poseidon",
        "",
        "babybear",
        24,
        "single",
        "8+21",
        Poseidon::new(&POSEIDON_BABYBEAR_24_PARAMS)
    );
}

fn add_poseidon2_cases(cases: &mut Vec<Case>) {
    add_case!(
        cases,
        "poseidon2",
        "",
        "bls12_381",
        3,
        "single",
        "8+56",
        Poseidon2::new(&POSEIDON2_BLS12_381_3_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "bls12_381",
        4,
        "single",
        "8+56",
        Poseidon2::new(&POSEIDON2_BLS12_381_4_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "bn254",
        3,
        "single",
        "8+56",
        Poseidon2::new(&POSEIDON2_BN254_3_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "bn254",
        4,
        "single",
        "8+56",
        Poseidon2::new(&POSEIDON2_BN254_4_GNARK_CRYPTO_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "goldilocks",
        8,
        "single",
        "8+22",
        Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "goldilocks",
        12,
        "single",
        "8+22",
        Poseidon2::new(&POSEIDON2_GOLDILOCKS_12_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "mersenne31",
        16,
        "single",
        "8+14",
        Poseidon2::new(&POSEIDON2_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "mersenne31",
        24,
        "single",
        "8+22",
        Poseidon2::new(&POSEIDON2_MERSENNE31_24_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "koalabear",
        16,
        "single",
        "8+20",
        Poseidon2::new(&POSEIDON2_KOALABEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "koalabear",
        24,
        "single",
        "8+23",
        Poseidon2::new(&POSEIDON2_KOALABEAR_24_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "babybear",
        16,
        "single",
        "8+13",
        Poseidon2::new(&POSEIDON2_BABYBEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "poseidon2",
        "",
        "babybear",
        24,
        "single",
        "8+21",
        Poseidon2::new(&POSEIDON2_BABYBEAR_24_PARAMS)
    );
}

fn add_type2_cases(cases: &mut Vec<Case>) {
    add_rescue_cases(cases);

    add_case!(
        cases,
        "anemoi",
        "",
        "bls12_381",
        4,
        "current",
        "20",
        Anemoi::new(&ANEMOI_BLS12_381_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "bls12_381",
        4,
        "original",
        "14",
        Anemoi::new(&ANEMOI_BLS12_381_4_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "bn254",
        4,
        "current",
        "20",
        Anemoi::new(&ANEMOI_BN254_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "bn254",
        4,
        "original",
        "14",
        Anemoi::new(&ANEMOI_BN254_4_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "goldilocks",
        8,
        "single",
        "11",
        Anemoi::new(&ANEMOI_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "goldilocks",
        12,
        "single",
        "10",
        Anemoi::new(&ANEMOI_GOLDILOCKS_12_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "mersenne31",
        16,
        "single",
        "10",
        Anemoi::new(&ANEMOI_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "mersenne31",
        24,
        "single",
        "9",
        Anemoi::new(&ANEMOI_MERSENNE31_24_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "koalabear",
        16,
        "single",
        "10",
        Anemoi::new(&ANEMOI_KOALABEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "koalabear",
        24,
        "single",
        "9",
        Anemoi::new(&ANEMOI_KOALABEAR_24_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "babybear",
        16,
        "single",
        "9",
        Anemoi::new(&ANEMOI_BABYBEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "anemoi",
        "",
        "babybear",
        24,
        "single",
        "9",
        Anemoi::new(&ANEMOI_BABYBEAR_24_PARAMS)
    );

    add_griffin_cases(cases);

    add_case!(
        cases,
        "arion",
        "",
        "bls12_381",
        3,
        "single",
        "11",
        Arion::new(&ARION_BLS12_381_3_PARAMS)
    );
    add_case!(
        cases,
        "arion",
        "",
        "bls12_381",
        4,
        "single",
        "10",
        Arion::new(&ARION_BLS12_381_4_PARAMS)
    );
    add_case!(
        cases,
        "arion",
        "",
        "bn254",
        3,
        "single",
        "11",
        Arion::new(&ARION_BN254_3_PARAMS)
    );
    add_case!(
        cases,
        "arion",
        "",
        "bn254",
        4,
        "single",
        "10",
        Arion::new(&ARION_BN254_4_PARAMS)
    );

    add_case!(
        cases,
        "grendel",
        "native",
        "bls12_381",
        3,
        "current",
        "24",
        Grendel::new(&GRENDEL_BLS12_381_3_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "grendel",
        "native",
        "bls12_381",
        3,
        "original",
        "14",
        Grendel::new(&GRENDEL_BLS12_381_3_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "grendel",
        "native",
        "bls12_381",
        4,
        "current",
        "23",
        Grendel::new(&GRENDEL_BLS12_381_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "grendel",
        "native",
        "bls12_381",
        4,
        "original",
        "12",
        Grendel::new(&GRENDEL_BLS12_381_4_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "grendel",
        "native",
        "bn254",
        3,
        "current",
        "24",
        Grendel::new(&GRENDEL_BN254_3_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "grendel",
        "native",
        "bn254",
        3,
        "original",
        "14",
        Grendel::new(&GRENDEL_BN254_3_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "grendel",
        "native",
        "bn254",
        4,
        "current",
        "23",
        Grendel::new(&GRENDEL_BN254_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "grendel",
        "native",
        "bn254",
        4,
        "original",
        "12",
        Grendel::new(&GRENDEL_BN254_4_ORIGINAL_PARAMS)
    );

    add_case!(
        cases,
        "xhash",
        "XHash8",
        "goldilocks",
        12,
        "single",
        "6",
        xhash8_goldilocks()
    );
    add_case!(
        cases,
        "xhash",
        "XHash12",
        "goldilocks",
        12,
        "single",
        "6",
        xhash12_goldilocks()
    );
    add_case!(
        cases,
        "xhash",
        "XHash16",
        "mersenne31",
        24,
        "single",
        "6",
        xhash16_m31()
    );
    add_case!(
        cases,
        "xhash",
        "XHash24",
        "mersenne31",
        24,
        "single",
        "6",
        xhash24_m31()
    );
}

fn add_rescue_cases(cases: &mut Vec<Case>) {
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bls12_381",
        3,
        "current",
        "16",
        RescuePrime::new(&RESCUE_PRIME_BLS12_381_3_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bls12_381",
        3,
        "original",
        "14",
        RescuePrime::new(&RESCUE_PRIME_BLS12_381_3_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bls12_381",
        4,
        "current",
        "13",
        RescuePrime::new(&RESCUE_PRIME_BLS12_381_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bls12_381",
        4,
        "original",
        "11",
        RescuePrime::new(&RESCUE_PRIME_BLS12_381_4_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bn254",
        3,
        "current",
        "16",
        RescuePrime::new(&RESCUE_PRIME_BN254_3_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bn254",
        3,
        "original",
        "14",
        RescuePrime::new(&RESCUE_PRIME_BN254_3_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bn254",
        4,
        "current",
        "13",
        RescuePrime::new(&RESCUE_PRIME_BN254_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "bn254",
        4,
        "original",
        "11",
        RescuePrime::new(&RESCUE_PRIME_BN254_4_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "goldilocks",
        8,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "goldilocks",
        12,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_GOLDILOCKS_12_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "mersenne31",
        16,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "mersenne31",
        24,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_MERSENNE31_24_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "koalabear",
        16,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_KOALABEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "koalabear",
        24,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_KOALABEAR_24_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "babybear",
        16,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_BABYBEAR_16_PARAMS)
    );
    add_case!(
        cases,
        "rescue_prime",
        "",
        "babybear",
        24,
        "single",
        "8",
        RescuePrime::new(&RESCUE_PRIME_BABYBEAR_24_PARAMS)
    );
}

fn add_griffin_cases(cases: &mut Vec<Case>) {
    add_case!(
        cases,
        "griffin",
        "",
        "bls12_381",
        3,
        "current",
        "20",
        Griffin::new(&GRIFFIN_BLS12_381_3_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "bls12_381",
        3,
        "original",
        "14",
        Griffin::new(&GRIFFIN_BLS12_381_3_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "bls12_381",
        4,
        "current",
        "20",
        Griffin::new(&GRIFFIN_BLS12_381_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "bls12_381",
        4,
        "original",
        "11",
        Griffin::new(&GRIFFIN_BLS12_381_4_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "bn254",
        3,
        "current",
        "20",
        Griffin::new(&GRIFFIN_BN254_3_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "bn254",
        3,
        "original",
        "14",
        Griffin::new(&GRIFFIN_BN254_3_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "bn254",
        4,
        "current",
        "20",
        Griffin::new(&GRIFFIN_BN254_4_CURRENT_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "bn254",
        4,
        "original",
        "11",
        Griffin::new(&GRIFFIN_BN254_4_ORIGINAL_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "goldilocks",
        8,
        "single",
        "8",
        Griffin::new(&GRIFFIN_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "griffin",
        "",
        "goldilocks",
        12,
        "single",
        "8",
        Griffin::new(&GRIFFIN_GOLDILOCKS_12_PARAMS)
    );
}

fn add_type3_cases(cases: &mut Vec<Case>) {
    add_case!(
        cases,
        "reinforced_concrete",
        "",
        "bls12_381",
        3,
        "single",
        "7",
        ReinforcedConcrete::new(&REINFORCED_CONCRETE_BLS12_381_3_PARAMS)
    );
    add_case!(
        cases,
        "reinforced_concrete",
        "",
        "bn254",
        3,
        "single",
        "7",
        ReinforcedConcrete::new(&REINFORCED_CONCRETE_BN254_3_PARAMS)
    );
    add_skyscraper(
        cases,
        "bls12_381",
        Skyscraper::new(&SKYSCRAPER_BLS12_381_2_PARAMS),
    );
    add_skyscraper(cases, "bn254", Skyscraper::new(&SKYSCRAPER_BN254_2_PARAMS));

    add_case!(
        cases,
        "monolith",
        "",
        "goldilocks",
        8,
        "single",
        "6",
        Monolith64::new(&MONOLITH_GOLDILOCKS_8_PARAMS)
    );
    add_case!(
        cases,
        "monolith",
        "",
        "goldilocks",
        12,
        "single",
        "6",
        Monolith64::new(&MONOLITH_GOLDILOCKS_12_PARAMS)
    );
    add_case!(
        cases,
        "monolith",
        "",
        "mersenne31",
        16,
        "single",
        "6",
        Monolith31::new(&MONOLITH_MERSENNE31_16_PARAMS)
    );
    add_case!(
        cases,
        "monolith",
        "",
        "mersenne31",
        24,
        "single",
        "6",
        Monolith31::new(&MONOLITH_MERSENNE31_24_PARAMS)
    );

    add_case!(
        cases,
        "polocolo",
        "",
        "bls12_381",
        3,
        "single",
        "6",
        Polocolo::new(&POLOCOLO_BLS12_381_3_PARAMS)
    );
    add_case!(
        cases,
        "polocolo",
        "",
        "bls12_381",
        4,
        "single",
        "5",
        Polocolo::new(&POLOCOLO_BLS12_381_4_PARAMS)
    );
    add_case!(
        cases,
        "polocolo",
        "",
        "bn254",
        3,
        "single",
        "6",
        Polocolo::new(&POLOCOLO_BN254_3_PARAMS)
    );
    add_case!(
        cases,
        "polocolo",
        "",
        "bn254",
        4,
        "single",
        "5",
        Polocolo::new(&POLOCOLO_BN254_4_PARAMS)
    );
    add_case!(
        cases,
        "tip4p",
        "",
        "goldilocks",
        12,
        "single",
        "5",
        Tip4::new(&TIP4P_GOLDILOCKS_PARAMS)
    );
}

fn add_skyscraper<F: PrimeFieldMontgomery + 'static>(
    cases: &mut Vec<Case>,
    field: &'static str,
    permutation: Skyscraper<F, 2>,
) {
    add_permutation(
        cases,
        ("skyscraper", "", field, 4, "single", "18"),
        move |input| permutation.permutation(input),
    );
}

fn add_plain_cases(cases: &mut Vec<Case>) {
    cases.push(plain_case("sha256_compress", run_sha256));
    cases.push(plain_case("keccak_f1600", run_keccak));
    cases.push(plain_case("blake2b_compress", run_blake2b));
    cases.push(plain_case("blake3_compress", run_blake3));
}

fn plain_case(variant: &'static str, run: fn(usize) -> u128) -> Case {
    Case {
        construction: "plain",
        variant,
        field: "-",
        t: None,
        round_set: "single",
        rounds: "-",
        run: Box::new(run),
    }
}

fn run_sha256(iters: usize) -> u128 {
    use sha2::digest::generic_array::GenericArray;
    let mut state = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bytes: [u8; 64] = core::array::from_fn(|i| (i as u8).wrapping_add(1));
    let block = GenericArray::clone_from_slice(&bytes);
    let start = Instant::now();
    for _ in 0..iters {
        plain_hashes::sha256_compress(
            black_box(&mut state),
            black_box(core::slice::from_ref(&block)),
        );
        black_box(&state);
    }
    let elapsed = start.elapsed().as_nanos();
    black_box(state);
    elapsed
}

fn run_keccak(iters: usize) -> u128 {
    let mut state: [u64; 25] = core::array::from_fn(|i| i as u64);
    let start = Instant::now();
    for _ in 0..iters {
        plain_hashes::keccak_f1600(black_box(&mut state));
        black_box(&state);
    }
    let elapsed = start.elapsed().as_nanos();
    black_box(state);
    elapsed
}

fn run_blake2b(iters: usize) -> u128 {
    let mut state = [
        0x6a09e667f3bcc908,
        0xbb67ae8584caa73b,
        0x3c6ef372fe94f82b,
        0xa54ff53a5f1d36f1,
        0x510e527fade682d1,
        0x9b05688c2b3e6c1f,
        0x1f83d9abfb41bd6b,
        0x5be0cd19137e2179,
    ];
    let block: [u64; 16] = core::array::from_fn(|i| i as u64);
    let start = Instant::now();
    for _ in 0..iters {
        plain_hashes::blake2b_compress(black_box(&mut state), black_box(&block), 0, 0, 0, 0);
        black_box(&state);
    }
    let elapsed = start.elapsed().as_nanos();
    black_box(state);
    elapsed
}

fn run_blake3(iters: usize) -> u128 {
    let chaining_value = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let block: [u8; 64] = core::array::from_fn(|i| i as u8);
    let start = Instant::now();
    for _ in 0..iters {
        black_box(plain_hashes::blake3_compress(
            black_box(&chaining_value),
            black_box(&block),
            64,
            0,
            0,
        ));
    }
    start.elapsed().as_nanos()
}

fn print_preamble(options: &Options) {
    preamble(
        "benchmark_mode",
        match options.preview_iters {
            Some(iters) => format!("preview; fixed_iters={iters}; not_for_publication"),
            None => "publication; auto_calibrated".to_owned(),
        },
    );
    preamble("cpu_model", cpu_model());
    preamble("rustc", command_output("rustc", &["--version"]));
    preamble("git_commit", command_output("git", &["rev-parse", "HEAD"]));
    let dirty = Command::new("git")
        .args(["status", "--porcelain", "--untracked-files=normal"])
        .output()
        .map(|output| {
            if output.stdout.is_empty() {
                "clean"
            } else {
                "dirty"
            }
            .to_owned()
        })
        .unwrap_or_else(|error| format!("unavailable ({error})"));
    preamble("git_state", dirty);
    preamble(
        "RUSTFLAGS",
        option_env!("RUSTFLAGS").unwrap_or("unset").to_owned(),
    );
    preamble(
        "release_profile",
        "lto=fat, codegen-units=1, panic=abort".to_owned(),
    );
    preamble("pinned_cpus", status_value("Cpus_allowed_list"));
    preamble("governor", governor());
    preamble("turbo", turbo_state());
    preamble("smt", read_trimmed("/sys/devices/system/cpu/smt/control"));
    preamble("command", env::args().collect::<Vec<_>>().join(" "));
}

fn preamble(key: &str, value: String) {
    println!("# {key}: {}", value.replace(['\n', '\r'], " "));
}

fn cpu_model() -> String {
    fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|contents| {
            contents.lines().find_map(|line| {
                let (key, value) = line.split_once(':')?;
                (key.trim() == "model name").then(|| value.trim().to_owned())
            })
        })
        .unwrap_or_else(|| "unavailable".to_owned())
}

fn command_output(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .map(|output| String::from_utf8_lossy(&output.stdout).trim().to_owned())
        .unwrap_or_else(|error| format!("unavailable ({error})"))
}

fn status_value(wanted: &str) -> String {
    fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|contents| {
            contents.lines().find_map(|line| {
                let (key, value) = line.split_once(':')?;
                (key == wanted).then(|| value.trim().to_owned())
            })
        })
        .unwrap_or_else(|| "unavailable".to_owned())
}

fn governor() -> String {
    let allowed = status_value("Cpus_allowed_list");
    let first_cpu = allowed
        .split(',')
        .next()
        .and_then(|part| part.split('-').next())
        .unwrap_or("0");
    read_trimmed(&format!(
        "/sys/devices/system/cpu/cpu{first_cpu}/cpufreq/scaling_governor"
    ))
}

fn turbo_state() -> String {
    for path in [
        "/sys/devices/system/cpu/intel_pstate/no_turbo",
        "/sys/devices/system/cpu/cpufreq/boost",
    ] {
        if let Ok(value) = fs::read_to_string(path) {
            return format!(
                "{}={}",
                path.rsplit('/').next().unwrap_or(path),
                value.trim()
            );
        }
    }
    "unavailable".to_owned()
}

fn read_trimmed(path: &str) -> String {
    fs::read_to_string(path)
        .map(|value| value.trim().to_owned())
        .unwrap_or_else(|error| format!("unavailable ({error})"))
}
