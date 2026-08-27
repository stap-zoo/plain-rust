pub mod babybear;
pub mod bls12_381;
pub mod bn254;
pub mod felt252;
pub mod goldilocks;
pub mod koalabear;
pub mod mersenne31;
mod montgomery_31;
mod montgomery_4;

use num_bigint::BigUint;
use num_traits::Zero;

pub trait FieldElement: Clone + Default + PartialEq + Eq + std::fmt::Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn from_u64(val: u64) -> Self;

    fn add_assign(&mut self, other: &Self);
    fn sub_assign(&mut self, other: &Self);
    fn mul_assign(&mut self, other: &Self);

    fn square(&mut self) {
        let tmp = self.clone();
        self.mul_assign(&tmp);
    }

    fn double(&mut self) {
        let tmp = self.clone();
        self.add_assign(&tmp);
    }

    fn negate(&self) -> Self {
        let mut out = Self::zero();
        out.sub_assign(self);
        out
    }

    fn pow_u64(&self, mut exp: u64) -> Self {
        let mut base = self.clone();
        let mut result = Self::one();
        while exp > 0 {
            if exp & 1 == 1 {
                result.mul_assign(&base);
            }
            exp >>= 1;
            if exp > 0 {
                base.square();
            }
        }
        result
    }

    /// Exponentiation where `exp_words_le` is represented as little-endian 64-bit limbs.
    fn pow_words_le(&self, exp_words_le: &[u64]) -> Self {
        let mut result = Self::one();
        let mut started = false;

        for &word in exp_words_le.iter().rev() {
            let mut mask = 1u64 << 63;
            while mask != 0 {
                if started {
                    result.square();
                }

                if (word & mask) != 0 {
                    if started {
                        result.mul_assign(self);
                    } else {
                        result = self.clone();
                        started = true;
                    }
                }

                mask >>= 1;
            }
        }

        result
    }
}

pub trait PrimeField: FieldElement {
    fn modulus() -> BigUint;
    fn from_biguint(value: &BigUint) -> Self;
    fn generator() -> BigUint;
}

pub trait PrimeFieldExt: PrimeField {
    fn to_biguint(&self) -> BigUint;
}

pub trait PrimeFieldWords: PrimeFieldExt {
    fn to_words_le(&self) -> [u64; 4];

    fn from_words_le(words: [u64; 4]) -> Self {
        Self::from_biguint(&biguint_from_limbs_le(&words))
    }
}

pub(crate) fn biguint_from_limbs_le(limbs: &[u64]) -> BigUint {
    let mut value = BigUint::zero();
    for limb in limbs.iter().rev() {
        value <<= 64;
        value += BigUint::from(*limb);
    }
    value
}

pub(crate) fn biguint_to_limbs_le_4(value: &BigUint) -> [u64; 4] {
    let mut out = [0u64; 4];
    let limbs = value.to_u64_digits();
    for (i, limb) in limbs.iter().enumerate().take(4) {
        out[i] = *limb;
    }
    out
}
