use super::{FieldElement, PrimeField, PrimeFieldExt, PrimeFieldMontgomery, PrimeFieldWords};
use ark_bls12_381::Fr as ArkBls12_381;
use ark_ff::{BigInt, BigInteger, PrimeField as ArkPrimeField};
use num_bigint::BigUint;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Bls12_381(pub(crate) ArkBls12_381);

impl Bls12_381 {
    pub fn from_hex(s: &str) -> Option<Self> {
        let s = s.strip_prefix("0x").unwrap_or(s);
        let s = if s.is_empty() { "0" } else { s };
        BigUint::parse_bytes(s.as_bytes(), 16).map(|v| Self(ArkBls12_381::from(v)))
    }
}

impl FieldElement for Bls12_381 {
    fn zero() -> Self {
        Self(ArkBls12_381::from(0u64))
    }

    fn one() -> Self {
        Self(ArkBls12_381::from(1u64))
    }

    fn from_u64(val: u64) -> Self {
        Self(ArkBls12_381::from(val))
    }

    fn add_assign(&mut self, other: &Self) {
        self.0 += other.0;
    }

    fn sub_assign(&mut self, other: &Self) {
        self.0 -= other.0;
    }

    fn mul_assign(&mut self, other: &Self) {
        self.0 *= other.0;
    }
}

impl PrimeField for Bls12_381 {
    fn modulus() -> BigUint {
        <ArkBls12_381 as ArkPrimeField>::MODULUS.into()
    }

    fn from_biguint(value: &BigUint) -> Self {
        Self(ArkBls12_381::from(value.clone()))
    }

    fn generator() -> BigUint {
        BigUint::from(7u32)
    }
}

impl PrimeFieldExt for Bls12_381 {
    fn to_biguint(&self) -> BigUint {
        self.0.into()
    }
}

impl PrimeFieldWords for Bls12_381 {
    fn to_words_le(&self) -> [u64; 4] {
        let limbs = self.0.into_bigint();
        let src = limbs.as_ref();
        let mut out = [0u64; 4];
        for (i, limb) in src.iter().enumerate().take(4) {
            out[i] = *limb;
        }
        out
    }

    fn from_words_le(words: [u64; 4]) -> Self {
        let mut bytes = [0u8; 32];
        for (chunk, word) in bytes.chunks_exact_mut(8).zip(words.iter()) {
            chunk.copy_from_slice(&word.to_le_bytes());
        }
        Self(ArkBls12_381::from_le_bytes_mod_order(&bytes))
    }
}

impl PrimeFieldMontgomery for Bls12_381 {
    fn into_montgomery_raw(self) -> Self {
        Self(ArkBls12_381::new_unchecked(self.0.into_bigint()))
    }

    fn from_montgomery_raw(self) -> Self {
        Self(ArkBls12_381::new(self.0 .0))
    }

    fn raw_words(&self) -> [u64; 4] {
        (self.0 .0).0
    }

    fn from_reduced_raw_words(words: [u64; 4]) -> Self {
        let mut value = BigInt(words);
        let modulus = <ArkBls12_381 as ArkPrimeField>::MODULUS;
        while value >= modulus {
            value.sub_with_borrow(&modulus);
        }
        Self(ArkBls12_381::new_unchecked(value))
    }
}
