use crate::fields::{FieldElement, PrimeField};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Signed, Zero};
use sha3::digest::XofReader;

/// Build a Cauchy MDS matrix of size t×t.
///
/// A Cauchy matrix C[i][j] = 1 / (x_i + y_j) where all x_i and y_j are
/// distinct and x_i + y_j ≠ 0.  Every square Cauchy matrix over a field
/// with distinct x_i, y_j is guaranteed to be MDS.
///
/// We use x_i = i+1, y_j = j+1 (both 1-indexed), so C[i][j] = 1/(i+j+2).
/// This is a standard choice also used in Poseidon parameter generation.
pub fn build_cauchy_mds<F: crate::fields::PrimeField>(t: usize) -> Vec<Vec<F>> {
    let modulus = F::modulus();
    let mut m = vec![vec![F::zero(); t]; t];
    for i in 0..t {
        for j in 0..t {
            // C[i][j] = 1 / (i + j + 2)    (i, j are 0-indexed, so +2)
            let denom = BigUint::from((i + j + 2) as u64);
            let inv = modinv(&denom, &modulus);
            m[i][j] = F::from_biguint(&inv);
        }
    }
    m
}

/// Build the transposed Vandermonde-echelon MDS matrix used by Grendel and RPO.
///
/// Starting from `V[i][j] = g^(i*j)` for a `t x 2t` Vandermonde matrix, this
/// reduces the left half to the identity, obtaining `(I | M^T)`, and returns
/// the transpose of the right half.
pub fn build_vandermonde_mds_transposed<F: PrimeField>(
    t: usize,
    generator: &BigUint,
) -> Vec<Vec<F>> {
    assert!(t > 0, "Vandermonde MDS width must be positive");

    let g = F::from_biguint(generator);
    let mut augmented = vec![vec![F::zero(); 2 * t]; t];
    for (i, row) in augmented.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = g.pow_u64((i * j) as u64);
        }
    }

    let modulus_minus_two = F::modulus() - BigUint::from(2u64);
    let inverse_exp = modulus_minus_two.to_u64_digits();

    for col in 0..t {
        let pivot = (col..t)
            .find(|&row| augmented[row][col] != F::zero())
            .expect("Vandermonde left half must be invertible");
        augmented.swap(col, pivot);

        let pivot_inv = augmented[col][col].pow_words_le(&inverse_exp);
        for entry in &mut augmented[col] {
            entry.mul_assign(&pivot_inv);
        }

        for row in 0..t {
            if row == col {
                continue;
            }
            let factor = augmented[row][col].clone();
            if factor == F::zero() {
                continue;
            }
            for j in 0..2 * t {
                let mut term = augmented[col][j].clone();
                term.mul_assign(&factor);
                augmented[row][j].sub_assign(&term);
            }
        }
    }

    let mut result = vec![vec![F::zero(); t]; t];
    for i in 0..t {
        for j in 0..t {
            result[i][j] = augmented[j][t + i].clone();
        }
    }
    result
}

pub(crate) fn modinv(value: &BigUint, modulus: &BigUint) -> BigUint {
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
        panic!("value is not invertible modulo modulus");
    }

    if t.is_negative() {
        t += BigInt::from(modulus.clone());
    }

    t.try_into().expect("modular inverse must be non-negative")
}

pub(crate) fn pow_biguint<F: FieldElement>(base: &F, exp: &BigUint) -> F {
    let mut result = F::one();
    let mut base_power = base.clone();
    let mut e = exp.clone();

    while !e.is_zero() {
        if (&e & BigUint::one()) == BigUint::one() {
            result.mul_assign(&base_power);
        }
        e >>= 1;
        if !e.is_zero() {
            base_power.square();
        }
    }

    result
}

/// Sample a uniformly random field element from a SHAKE128 XOF stream.
///
/// Uses masked-byte rejection sampling: reads `ceil(modulus_bits/8)` bytes,
/// masks the top bits of the last byte, and rejects values ≥ modulus.
/// This is the standard method used by zkfriendlyhashzoo.
pub fn read_field_from_shake<F: PrimeField>(reader: &mut dyn XofReader) -> F {
    let modulus = F::modulus();
    let bits = modulus.bits() as usize;
    let byte_len = (bits + 7) / 8;
    let mod_bits = bits % 8;
    let mask = if mod_bits == 0 {
        0xFFu8
    } else {
        (1u8 << mod_bits) - 1
    };
    let last = byte_len.saturating_sub(1);
    let mut buf = vec![0u8; byte_len];

    loop {
        reader.read(&mut buf);
        if mod_bits != 0 {
            buf[last] &= mask;
        }
        let val = BigUint::from_bytes_le(&buf);
        if val < modulus {
            return F::from_biguint(&val);
        }
    }
}

/// Read one big-endian SHAKE chunk and reduce it modulo the field modulus.
///
/// Grendel uses `w = ceil(modulus_bits / 8) + 1`, but the width remains an
/// argument so this helper also covers other modulo-sampling specifications.
pub fn read_field_from_shake256_mod<F: PrimeField>(reader: &mut dyn XofReader, width: usize) -> F {
    assert!(width > 0, "SHAKE field-element chunk must be non-empty");
    let mut bytes = vec![0u8; width];
    reader.read(&mut bytes);
    F::from_biguint(&(BigUint::from_bytes_be(&bytes) % F::modulus()))
}
