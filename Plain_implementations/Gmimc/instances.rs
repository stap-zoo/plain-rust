use super::gmimc::GmimcParams;
use crate::fields::babybear::BabyBear;
use crate::fields::bls12_381::Bls12_381;
use crate::fields::bn254::Bn254;
use crate::fields::goldilocks::Goldilocks;
use crate::fields::koalabear::KoalaBear;
use crate::fields::mersenne31::Mersenne31;
use crate::fields::PrimeField;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use std::sync::Arc;

/// Derive the reference constants with
/// `SHAKE256("GMiMC(p,t,R)aff")`, sampling one little-endian field element
/// from `ceil(bitlen(p)/8) + 1` bytes and reducing modulo `p`.
fn derive_params<F: PrimeField>(t: usize, alpha: u64, rounds: usize) -> Arc<GmimcParams<F>> {
    let modulus = F::modulus();
    let seed = format!("GMiMC({modulus},{t},{rounds})aff");
    let mut shake = Shake256::default();
    shake.update(seed.as_bytes());
    let mut reader = shake.finalize_xof();

    let draw_bytes = (modulus.bits() as usize).div_ceil(8) + 1;
    let mut round_constants = Vec::with_capacity(rounds);
    for _ in 0..rounds {
        let mut bytes = vec![0u8; draw_bytes];
        reader.read(&mut bytes);
        round_constants.push(F::from_biguint(&BigUint::from_bytes_le(&bytes)));
    }

    Arc::new(GmimcParams::new(t, alpha, rounds, &round_constants))
}

lazy_static! {
    // log2(p) ~= 64: alpha 7 from the Plonky3 benchmark plan.
    pub static ref GMIMC_GOLDILOCKS_8_PARAMS: Arc<GmimcParams<Goldilocks>> =
        derive_params(8, 7, 68);
    pub static ref GMIMC_GOLDILOCKS_12_PARAMS: Arc<GmimcParams<Goldilocks>> =
        derive_params(12, 7, 93);

    // log2(p) ~= 32: field-specific alpha from the Plonky3 benchmark plan.
    pub static ref GMIMC_MERSENNE31_16_PARAMS: Arc<GmimcParams<Mersenne31>> =
        derive_params(16, 5, 158);
    pub static ref GMIMC_MERSENNE31_24_PARAMS: Arc<GmimcParams<Mersenne31>> =
        derive_params(24, 5, 335);
    pub static ref GMIMC_BABYBEAR_16_PARAMS: Arc<GmimcParams<BabyBear>> =
        derive_params(16, 7, 158);
    pub static ref GMIMC_BABYBEAR_24_PARAMS: Arc<GmimcParams<BabyBear>> =
        derive_params(24, 7, 335);
    pub static ref GMIMC_KOALABEAR_16_PARAMS: Arc<GmimcParams<KoalaBear>> =
        derive_params(16, 3, 158);
    pub static ref GMIMC_KOALABEAR_24_PARAMS: Arc<GmimcParams<KoalaBear>> =
        derive_params(24, 3, 335);

    // log2(p) ~= 256: only the requested ~768-bit states; alpha 5 and R=228
    // are the BN254/BLS12-381 choices registered in the Gnark benchmark plan.
    pub static ref GMIMC_BN254_3_PARAMS: Arc<GmimcParams<Bn254>> =
        derive_params(3, 5, 228);
    pub static ref GMIMC_BLS12_381_3_PARAMS: Arc<GmimcParams<Bls12_381>> =
        derive_params(3, 5, 228);
    pub static ref GMIMC_BN254_4_PARAMS: Arc<GmimcParams<Bn254>> =
        derive_params(4, 5, 231);
    pub static ref GMIMC_BLS12_381_4_PARAMS: Arc<GmimcParams<Bls12_381>> =
        derive_params(4, 5, 231);
}
