use super::arion::ArionParams;
use crate::fields::bls12_381::Bls12_381;
use crate::fields::bn254::Bn254;
use crate::fields::PrimeField;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use std::sync::Arc;

fn sample_mod_le<F: PrimeField>(reader: &mut dyn XofReader) -> F {
    let modulus = F::modulus();
    let width = (modulus.bits() as usize).div_ceil(8) + 1;
    let mut bytes = vec![0u8; width];
    reader.read(&mut bytes);
    F::from_biguint(&(BigUint::from_bytes_le(&bytes) % modulus))
}

fn shake_reader(seed: &str) -> impl XofReader {
    let mut shake = Shake256::default();
    shake.update(seed.as_bytes());
    shake.finalize_xof()
}

fn generate_params<F: PrimeField>(t: usize, rounds: usize) -> Arc<ArionParams<F>> {
    let modulus = F::modulus();
    let modulus_minus_one = &modulus - BigUint::from(1u64);
    let d = 257u64;
    let d_inv = compute_d_inv(d, &modulus_minus_one);
    let base = format!("Arion({modulus},{t},{rounds})");

    let mut aff_reader = shake_reader(&(base.clone() + "aff"));
    let round_constants: Vec<Vec<F>> = (0..rounds)
        .map(|_| (0..t).map(|_| sample_mod_le(&mut aff_reader)).collect())
        .collect();

    let mut h_reader = shake_reader(&(base.clone() + "h"));
    let h_coeffs: Vec<Vec<F>> = (0..rounds)
        .map(|_| (0..t - 1).map(|_| sample_mod_le(&mut h_reader)).collect())
        .collect();

    let needed = rounds * (t - 1);
    let mut g_reader = shake_reader(&(base + "g"));
    let mut accepted = Vec::with_capacity(needed);
    for _ in 0..4 * needed {
        let a = sample_mod_le(&mut g_reader);
        let b = sample_mod_le(&mut g_reader);
        if is_irreducible(&a, &b, &modulus, &modulus_minus_one) {
            accepted.push([a, b]);
            if accepted.len() == needed {
                break;
            }
        }
    }
    assert_eq!(
        accepted.len(),
        needed,
        "Arion coefficient candidate pool exhausted"
    );
    let g_coeffs = accepted
        .chunks_exact(t - 1)
        .map(|row| row.to_vec())
        .collect::<Vec<_>>();

    Arc::new(ArionParams::new(
        t,
        d,
        d_inv,
        rounds,
        &g_coeffs,
        &h_coeffs,
        &round_constants,
    ))
}

fn compute_d_inv(d: u64, modulus_minus_one: &num_bigint::BigUint) -> [u64; 4] {
    let d_big = num_bigint::BigUint::from(d);
    let inv = modinv_biguint(&d_big, modulus_minus_one);
    let mut out = [0u64; 4];
    let limbs = inv.to_u64_digits();
    for (i, limb) in limbs.iter().enumerate().take(4) {
        out[i] = *limb;
    }
    out
}

fn modinv_biguint(
    value: &num_bigint::BigUint,
    modulus: &num_bigint::BigUint,
) -> num_bigint::BigUint {
    use num_bigint::BigInt;
    use num_traits::{One, Signed, Zero};
    let mut t = BigInt::zero();
    let mut new_t = BigInt::one();
    let mut r = BigInt::from(modulus.clone());
    let mut new_r = BigInt::from(value.clone());
    while !new_r.is_zero() {
        let quotient = &r / &new_r;
        let next_t = &t - &quotient * &new_t;
        t = new_t;
        new_t = next_t;
        let next_r = &r - &quotient * &new_r;
        r = new_r;
        new_r = next_r;
    }
    if r != BigInt::one() {
        panic!("not invertible");
    }
    if t.is_negative() {
        t += BigInt::from(modulus.clone());
    }
    t.try_into().expect("non-negative")
}

/// Check if x^2 + a1*x + a2 is irreducible over F_p.
/// True iff discriminant D = a1^2 - 4*a2 is a quadratic non-residue.
fn is_irreducible<F: PrimeField>(
    a1: &F,
    a2: &F,
    _modulus: &num_bigint::BigUint,
    modulus_minus_one: &num_bigint::BigUint,
) -> bool {
    use num_bigint::BigUint;
    // Discriminant D = a1^2 - 4*a2
    let mut a1_sq = a1.clone();
    a1_sq.square();
    let mut four_a2 = a2.clone();
    four_a2.mul_assign(&F::from_u64(4));
    let mut d = a1_sq;
    d.sub_assign(&four_a2);

    // D == 0 => double root => reducible
    if d == F::zero() {
        return false;
    }

    // Euler's criterion: Legendre(D) = D^{(p-1)/2} mod p.
    // Irreducible iff Legendre(D) = -1 (quadratic non-residue),
    // i.e. result == p - 1.
    let half = modulus_minus_one / BigUint::from(2u64);
    let half_limbs = half.to_u64_digits();
    let d_legendre = d.pow_words_le(&half_limbs[..]);

    let p_minus_one = F::from_biguint(modulus_minus_one);
    d_legendre == p_minus_one
}

// ---------------------------------------------------------------------------
// Pre-computed lazy_static instances
// ---------------------------------------------------------------------------

lazy_static! {
    // -- BN254, t=3 --
    pub static ref ARION_BN254_3_PARAMS: Arc<ArionParams<Bn254>> =
        generate_params::<Bn254>(3, 11);

    // -- BLS12-381, t=3 --
    pub static ref ARION_BLS12_381_3_PARAMS: Arc<ArionParams<Bls12_381>> =
        generate_params::<Bls12_381>(3, 11);
    pub static ref ARION_BN254_4_PARAMS: Arc<ArionParams<Bn254>> =
        generate_params::<Bn254>(4, 10);
    pub static ref ARION_BLS12_381_4_PARAMS: Arc<ArionParams<Bls12_381>> =
        generate_params::<Bls12_381>(4, 10);
}
