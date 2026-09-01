use super::grendel::GrendelParams;
use crate::fields::babybear::BabyBear;
use crate::fields::bls12_381::Bls12_381;
use crate::fields::bn254::Bn254;
use crate::fields::goldilocks::Goldilocks;
use crate::fields::koalabear::KoalaBear;
use crate::fields::mersenne31::Mersenne31;
use crate::fields::PrimeField;
use crate::utils::{build_vandermonde_mds_transposed, read_field_from_shake256_mod};
use lazy_static::lazy_static;
use sha3::digest::{ExtendableOutput, Update};
use sha3::Shake256;
use std::sync::Arc;

const KAPPA: usize = 128;

fn generate_params<F: PrimeField>(t: usize, alpha: u64, rounds: usize) -> Arc<GrendelParams<F>> {
    let modulus = F::modulus();
    let mds = build_vandermonde_mds_transposed::<F>(t, &F::generator());

    let seed = format!("grendel-{modulus}-{t}-{KAPPA}");
    let mut shake = Shake256::default();
    shake.update(seed.as_bytes());
    let mut reader = shake.finalize_xof();
    let width = ((modulus.bits() as usize + 7) / 8) + 1;
    let round_constants = (0..rounds)
        .map(|_| {
            (0..t)
                .map(|_| read_field_from_shake256_mod::<F>(&mut reader, width))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    Arc::new(GrendelParams::new(t, alpha, rounds, &mds, &round_constants))
}

lazy_static! {
    // The BLS12-381 and BN254 benchmark profiles use the same native exponent
    // for both widths. Current and original profiles differ only in rounds.
    pub static ref GRENDEL_BN254_3_CURRENT_PARAMS: Arc<GrendelParams<Bn254>> =
        generate_params::<Bn254>(3, 5, 24);
    pub static ref GRENDEL_BN254_3_ORIGINAL_PARAMS: Arc<GrendelParams<Bn254>> =
        generate_params::<Bn254>(3, 5, 14);
    pub static ref GRENDEL_BLS12_381_3_CURRENT_PARAMS: Arc<GrendelParams<Bls12_381>> =
        generate_params::<Bls12_381>(3, 5, 24);
    pub static ref GRENDEL_BLS12_381_3_ORIGINAL_PARAMS: Arc<GrendelParams<Bls12_381>> =
        generate_params::<Bls12_381>(3, 5, 14);
    pub static ref GRENDEL_BN254_4_CURRENT_PARAMS: Arc<GrendelParams<Bn254>> =
        generate_params::<Bn254>(4, 5, 23);
    pub static ref GRENDEL_BN254_4_ORIGINAL_PARAMS: Arc<GrendelParams<Bn254>> =
        generate_params::<Bn254>(4, 5, 12);
    pub static ref GRENDEL_BLS12_381_4_CURRENT_PARAMS: Arc<GrendelParams<Bls12_381>> =
        generate_params::<Bls12_381>(4, 5, 23);
    pub static ref GRENDEL_BLS12_381_4_ORIGINAL_PARAMS: Arc<GrendelParams<Bls12_381>> =
        generate_params::<Bls12_381>(4, 5, 12);
    pub static ref GRENDEL_GOLDILOCKS_8_PARAMS: Arc<GrendelParams<Goldilocks>> =
        generate_params::<Goldilocks>(8, 7, 15);
    pub static ref GRENDEL_GOLDILOCKS_12_PARAMS: Arc<GrendelParams<Goldilocks>> =
        generate_params::<Goldilocks>(12, 7, 11);
    pub static ref GRENDEL_MERSENNE31_16_PARAMS: Arc<GrendelParams<Mersenne31>> =
        generate_params::<Mersenne31>(16, 2, 12);
    pub static ref GRENDEL_MERSENNE31_24_PARAMS: Arc<GrendelParams<Mersenne31>> =
        generate_params::<Mersenne31>(24, 2, 12);
    pub static ref GRENDEL_BABYBEAR_16_PARAMS: Arc<GrendelParams<BabyBear>> =
        generate_params::<BabyBear>(16, 7, 14);
    pub static ref GRENDEL_BABYBEAR_24_PARAMS: Arc<GrendelParams<BabyBear>> =
        generate_params::<BabyBear>(24, 7, 14);
    pub static ref GRENDEL_KOALABEAR_16_PARAMS: Arc<GrendelParams<KoalaBear>> =
        generate_params::<KoalaBear>(16, 3, 12);
    pub static ref GRENDEL_KOALABEAR_24_PARAMS: Arc<GrendelParams<KoalaBear>> =
        generate_params::<KoalaBear>(24, 3, 12);
}

// Preserve the former current-profile names for downstream users.
pub use GRENDEL_BLS12_381_3_CURRENT_PARAMS as GRENDEL_BLS12_381_3_PARAMS;
pub use GRENDEL_BN254_3_CURRENT_PARAMS as GRENDEL_BN254_3_PARAMS;
