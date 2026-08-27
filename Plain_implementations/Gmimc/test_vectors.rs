//! Permutation KATs generated through the shared Python reference used by the
//! Plonky3 and Gnark sister repositories. Existing sister-repository vectors
//! were checked byte-for-byte at every overlapping instance. Each instance is
//! tested on the counting state `1..=t` and the all-zero state.

use super::gmimc::{Gmimc, GmimcParams};
use super::instances::*;
use crate::fields::babybear::BabyBear;
use crate::fields::bls12_381::Bls12_381;
use crate::fields::bn254::Bn254;
use crate::fields::goldilocks::Goldilocks;
use crate::fields::koalabear::KoalaBear;
use crate::fields::mersenne31::Mersenne31;
use crate::fields::PrimeField;
use num_bigint::BigUint;
use std::sync::Arc;

fn field_hex<F: PrimeField>(value: &str) -> F {
    let digits = value.strip_prefix("0x").unwrap_or(value);
    F::from_biguint(&BigUint::parse_bytes(digits.as_bytes(), 16).expect("valid KAT hex"))
}

fn check_vectors<F: PrimeField>(
    params: &Arc<GmimcParams<F>>,
    counting_output: &[&str],
    zero_output: &[&str],
) {
    let permutation = Gmimc::new(params);
    let t = permutation.get_t();
    assert_eq!(counting_output.len(), t);
    assert_eq!(zero_output.len(), t);

    let counting: Vec<F> = (1..=t).map(|i| F::from_u64(i as u64)).collect();
    let zero = vec![F::zero(); t];
    let expected_counting: Vec<F> = counting_output.iter().map(|x| field_hex(x)).collect();
    let expected_zero: Vec<F> = zero_output.iter().map(|x| field_hex(x)).collect();

    assert_eq!(permutation.permutation(&counting), expected_counting);
    assert_eq!(permutation.permutation(&zero), expected_zero);

    // The accumulator path is separately checked against the literal round loop.
    assert_eq!(
        permutation.permutation(&counting),
        permutation.permutation_reference(&counting)
    );
    assert_eq!(
        permutation.permutation(&zero),
        permutation.permutation_reference(&zero)
    );
}

#[test]
fn goldilocks_t8_vectors() {
    check_vectors::<Goldilocks>(
        &GMIMC_GOLDILOCKS_8_PARAMS,
        &[
            "0xf8b039686bb1c39f",
            "0xa7dcaa94efd3cdf5",
            "0xc38e347f9284e2ce",
            "0xfef88c98160a1c12",
            "0xc6172632d61f6798",
            "0x657451de5923c38",
            "0x5bfb3b3c269d3493",
            "0x3107e102795729f7",
        ],
        &[
            "0xc1b9f742c14d822d",
            "0x27dbe2bb853f90a",
            "0x6fddecbf48b5316c",
            "0xf815b4b5b6dd682e",
            "0xf03bcfc69cb670b7",
            "0xcff38a799a8115a1",
            "0xa372af6e85df8fa7",
            "0xc45b6e2a4a4d8b5a",
        ],
    );
}

#[test]
fn goldilocks_t12_vectors() {
    check_vectors::<Goldilocks>(
        &GMIMC_GOLDILOCKS_12_PARAMS,
        &[
            "0xf00694f64deb6f75",
            "0x13d3e9b8d667a2",
            "0x9718cd4465252a40",
            "0x533861575b399ab",
            "0xbaf2a8430b6ec8e4",
            "0x2947e722c2c89566",
            "0x6341b79b14f443ac",
            "0xb2234b243d14a5ab",
            "0xea7c87dc7b9ffd22",
            "0x41adcebfe94fcce9",
            "0xa972460ce6444c53",
            "0x2364467351aedc1a",
        ],
        &[
            "0xfca62aec5ee4ef1b",
            "0x21326615a81b9dd3",
            "0x94e58a9989b65ba3",
            "0xfd355171f3e0bada",
            "0x986ec2aae17bb518",
            "0x6cf6b8d4e221e7dc",
            "0x7ff5e03ea9e6b090",
            "0x4c072c1172576d87",
            "0x4a3369c8b627c867",
            "0x53d3972cfc849669",
            "0xca869898415956e8",
            "0x62a617eb48e29ad8",
        ],
    );
}

#[test]
fn mersenne31_t16_vectors() {
    check_vectors::<Mersenne31>(
        &GMIMC_MERSENNE31_16_PARAMS,
        &[
            "0x46d446d0",
            "0x7fe1d715",
            "0x3a077885",
            "0x3cc56d0b",
            "0x608365f1",
            "0x56d94827",
            "0x25ab48c1",
            "0x295ab4f9",
            "0x3d380f61",
            "0x3d00fc96",
            "0x37c558b5",
            "0x549f4fa0",
            "0x40823fb5",
            "0x6fb4fe75",
            "0x3cbb765f",
            "0x56bfa171",
        ],
        &[
            "0x741bbee6",
            "0x30322fc2",
            "0x3c253f3f",
            "0x62ea7782",
            "0x65e16fe7",
            "0x5dec7bb9",
            "0x6fedeb91",
            "0x15e4aea3",
            "0x4a45f8e1",
            "0x2d3b9d1b",
            "0x6244108e",
            "0x4dd69b78",
            "0x6ca3f818",
            "0x48936759",
            "0x4c2f4f9a",
            "0x37a97f9d",
        ],
    );
}

#[test]
fn mersenne31_t24_vectors() {
    check_vectors::<Mersenne31>(
        &GMIMC_MERSENNE31_24_PARAMS,
        &[
            "0x7a519100",
            "0xee831e5",
            "0x3ff4d2a9",
            "0x3b5618d9",
            "0x26e5825d",
            "0x6257d850",
            "0x3fc5001b",
            "0x13d8c170",
            "0x42b03d96",
            "0x37a5f371",
            "0x1e29a93d",
            "0x1569399d",
            "0x2df2101c",
            "0x52009099",
            "0x315ed564",
            "0x595f9cbb",
            "0x4cd129fa",
            "0x6aaf370a",
            "0x5ce758a9",
            "0x17822e0d",
            "0x1484d5ea",
            "0x4d75823a",
            "0x79c98f97",
            "0x44a3bc0f",
        ],
        &[
            "0x68c0c71d",
            "0x72e0dbb1",
            "0x6ec76359",
            "0xc05abed",
            "0x326ea24b",
            "0x1645770b",
            "0x5d4c2906",
            "0x60cc40fd",
            "0x639722b9",
            "0x3e75abdf",
            "0x654072e4",
            "0x46c08b41",
            "0x794f2413",
            "0x20506cf5",
            "0x2e106ba2",
            "0x11031522",
            "0x44c3dc14",
            "0x5423a153",
            "0x4cf82c2d",
            "0x5172d0d1",
            "0x1a79cb3a",
            "0x29003cfd",
            "0x653af5e1",
            "0x4ccce3af",
        ],
    );
}

#[test]
fn babybear_t16_vectors() {
    check_vectors::<BabyBear>(
        &GMIMC_BABYBEAR_16_PARAMS,
        &[
            "0x77933e89",
            "0x64054df9",
            "0x2c223e16",
            "0x2333650e",
            "0xb2e43c9",
            "0x708be1e0",
            "0x39363f5f",
            "0x4ffea5b5",
            "0x149eecc0",
            "0x2bc9fa4b",
            "0x118a8403",
            "0x1e1dd369",
            "0x5aa20fc",
            "0xfb3ea77",
            "0x47e99035",
            "0x5a547c37",
        ],
        &[
            "0x75452dc2",
            "0x862ed11",
            "0x33fc0981",
            "0x64991c42",
            "0x34499066",
            "0x2fec36f1",
            "0x2e3c2243",
            "0x2e88969f",
            "0x4f3be34b",
            "0x24cc6d00",
            "0x6b628c09",
            "0x1ea05255",
            "0x21068bf1",
            "0x1e46ae1a",
            "0x73935b94",
            "0x5eb66245",
        ],
    );
}

#[test]
fn babybear_t24_vectors() {
    check_vectors::<BabyBear>(
        &GMIMC_BABYBEAR_24_PARAMS,
        &[
            "0x7598b001",
            "0x728aadcb",
            "0x6c55cf1e",
            "0x94a3d19",
            "0x2a09bab2",
            "0x40144cef",
            "0x73cd60c9",
            "0x292ec8a2",
            "0x7518bac",
            "0x39ee6f4f",
            "0x1a8b5957",
            "0x5fff4383",
            "0x178ddc95",
            "0x4dbd76e4",
            "0x1351e56",
            "0x4fa1e681",
            "0x5e582079",
            "0x417a9ae7",
            "0x2c6cfed8",
            "0x74ea9bc2",
            "0x6ab7473f",
            "0x3d96e49c",
            "0x6b1b5ee7",
            "0x15764ca5",
        ],
        &[
            "0x6a5fb765",
            "0x27f58f16",
            "0x1f94aa00",
            "0x4d87986b",
            "0x2a0a73e2",
            "0x7574fdf1",
            "0x62904b3b",
            "0x28c49b2",
            "0x652a77f5",
            "0x5ffc8a0c",
            "0x1bf40d43",
            "0x6f1097fb",
            "0x3d0902c1",
            "0x3b7a0a45",
            "0x5eb8aca5",
            "0x710ce9de",
            "0x51f6ec44",
            "0x133de814",
            "0x3734d7e8",
            "0x2892137b",
            "0x57edbfdb",
            "0x2d433617",
            "0x5246000f",
            "0x61d178f5",
        ],
    );
}

#[test]
fn koalabear_t16_vectors() {
    check_vectors::<KoalaBear>(
        &GMIMC_KOALABEAR_16_PARAMS,
        &[
            "0x1ec90e12",
            "0x4eba3889",
            "0x666a755e",
            "0x2841e7c0",
            "0x661e613a",
            "0x1e5fbd50",
            "0x1e8aca7f",
            "0x76b78086",
            "0x7b2ef031",
            "0x5ce6e572",
            "0x75052482",
            "0x73b97118",
            "0x3409b646",
            "0x1b074e9",
            "0x1c7cbedd",
            "0x8908713",
        ],
        &[
            "0x2651f6b8",
            "0x37980c8a",
            "0x288edd64",
            "0x35805458",
            "0x5dec6997",
            "0x5a95154a",
            "0x23f7c3cf",
            "0x1c22caf",
            "0x41e90d41",
            "0x3f22be9c",
            "0x25cc61",
            "0x72764a54",
            "0x3465a33a",
            "0x215eac35",
            "0x8e5329f",
            "0x7bf1377f",
        ],
    );
}

#[test]
fn koalabear_t24_vectors() {
    check_vectors::<KoalaBear>(
        &GMIMC_KOALABEAR_24_PARAMS,
        &[
            "0xbfdb251",
            "0x67231621",
            "0x4a2892b7",
            "0x2f27a021",
            "0x41f85f68",
            "0x14e2b9ff",
            "0x37b6c813",
            "0x259b697d",
            "0x47359a8e",
            "0x40af0e05",
            "0x6b9d402",
            "0x2599346a",
            "0xc536df4",
            "0x66dda1c2",
            "0x305b72d0",
            "0x208f7adb",
            "0x392f2af8",
            "0x7cebfefd",
            "0x5896d8d0",
            "0x4faa2d65",
            "0xe043a9c",
            "0x82349ba",
            "0x1312162f",
            "0xede5f00",
        ],
        &[
            "0x156ce924",
            "0x681d154f",
            "0xa38de65",
            "0x2360b7e5",
            "0x1fe7dbc4",
            "0x65f2e019",
            "0x33ab637f",
            "0x50912dc",
            "0x78f6f597",
            "0x1c7fb47f",
            "0x3db6bf22",
            "0x13289ebd",
            "0x5adbb645",
            "0x3dd7a4f6",
            "0x511b5623",
            "0x24d5cf5c",
            "0x494f6c58",
            "0x60c4c868",
            "0x5b64e115",
            "0x27c0677",
            "0x375156f4",
            "0x3fb9c23a",
            "0x79d8e8c7",
            "0x6ca7b830",
        ],
    );
}

#[test]
fn bn254_t3_vectors() {
    check_vectors::<Bn254>(
        &GMIMC_BN254_3_PARAMS,
        &[
            "0x379ad8aac8ac2b87d5005b441c1956febf8e80830788afffc4215e5110457e3",
            "0x1d1076329a0d81f1539076dcb27e647566b964f1d05b31f60e1e0009632834f0",
            "0x25067cb46cb4072d2861d8c5b99651ec3b86219bbe9888d661fcf774cafade09",
        ],
        &[
            "0x15786830a9a1ed7ed344bd97a300d4a46961db6c49e6cbefb832b656c40bb331",
            "0x1176103d11445c7585bdef3d6fa19834d575c72b1f76cc4afb3e066216801fc0",
            "0xd3ad3d767359d5a883884093fbbb15353cae474a9525222e8073bf97c1d923c",
        ],
    );
}

#[test]
fn bls12_381_t3_vectors() {
    check_vectors::<Bls12_381>(
        &GMIMC_BLS12_381_3_PARAMS,
        &[
            "0x45fe6982a6c1dc21f9de5f7d0e9042efbab2e943efdae8a0b8483b116dc08b1b",
            "0x6496508f1d2609c60a2e481941438d047ccd2750a94d4c7e757b7caeea513208",
            "0x4373f5044782b7141dd9d4e55fadb7d50b6a3104681639940ca96b328654b965",
        ],
        &[
            "0x4883bfe5f65f1c9404bab42f08db5adae72c3dba87ca9bf690cf4f2be6f91298",
            "0x8a7d867b99143266794dea8572d3d976817896190c434d705b85d91ac6ccde4",
            "0x2213aba302d15d33d6e8202d75706297185d6813aa2e2e3da89925fcd120c315",
        ],
    );
}

#[test]
fn instance_shapes_match_the_round_table_and_benchmark_alphas() {
    let small = [
        (
            GMIMC_GOLDILOCKS_8_PARAMS.get_t(),
            GMIMC_GOLDILOCKS_8_PARAMS.get_alpha(),
            GMIMC_GOLDILOCKS_8_PARAMS.get_rounds(),
        ),
        (
            GMIMC_GOLDILOCKS_12_PARAMS.get_t(),
            GMIMC_GOLDILOCKS_12_PARAMS.get_alpha(),
            GMIMC_GOLDILOCKS_12_PARAMS.get_rounds(),
        ),
    ];
    assert_eq!(small, [(8, 7, 68), (12, 7, 93)]);

    assert_eq!(
        (
            GMIMC_MERSENNE31_16_PARAMS.get_alpha(),
            GMIMC_MERSENNE31_16_PARAMS.get_rounds()
        ),
        (5, 158)
    );
    assert_eq!(
        (
            GMIMC_MERSENNE31_24_PARAMS.get_alpha(),
            GMIMC_MERSENNE31_24_PARAMS.get_rounds()
        ),
        (5, 335)
    );
    assert_eq!(
        (
            GMIMC_BABYBEAR_16_PARAMS.get_alpha(),
            GMIMC_BABYBEAR_16_PARAMS.get_rounds()
        ),
        (7, 158)
    );
    assert_eq!(
        (
            GMIMC_BABYBEAR_24_PARAMS.get_alpha(),
            GMIMC_BABYBEAR_24_PARAMS.get_rounds()
        ),
        (7, 335)
    );
    assert_eq!(
        (
            GMIMC_KOALABEAR_16_PARAMS.get_alpha(),
            GMIMC_KOALABEAR_16_PARAMS.get_rounds()
        ),
        (3, 158)
    );
    assert_eq!(
        (
            GMIMC_KOALABEAR_24_PARAMS.get_alpha(),
            GMIMC_KOALABEAR_24_PARAMS.get_rounds()
        ),
        (3, 335)
    );
    assert_eq!(
        (
            GMIMC_BN254_3_PARAMS.get_t(),
            GMIMC_BN254_3_PARAMS.get_alpha(),
            GMIMC_BN254_3_PARAMS.get_rounds()
        ),
        (3, 5, 228)
    );
    assert_eq!(
        (
            GMIMC_BLS12_381_3_PARAMS.get_t(),
            GMIMC_BLS12_381_3_PARAMS.get_alpha(),
            GMIMC_BLS12_381_3_PARAMS.get_rounds()
        ),
        (3, 5, 228)
    );
}
