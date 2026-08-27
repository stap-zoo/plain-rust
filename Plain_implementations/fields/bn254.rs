use super::{FieldElement, PrimeField, PrimeFieldExt, PrimeFieldWords};
use ark_bn254::Fr as ArkBn254;
use ark_ff::PrimeField as ArkPrimeField;
use num_bigint::BigUint;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Bn254(pub(crate) ArkBn254);

impl Bn254 {
    pub fn from_hex(s: &str) -> Option<Self> {
        let s = s.strip_prefix("0x").unwrap_or(s);
        let s = if s.is_empty() { "0" } else { s };
        BigUint::parse_bytes(s.as_bytes(), 16).map(|v| Self(ArkBn254::from(v)))
    }
}

impl FieldElement for Bn254 {
    fn zero() -> Self {
        Self(ArkBn254::from(0u64))
    }

    fn one() -> Self {
        Self(ArkBn254::from(1u64))
    }

    fn from_u64(val: u64) -> Self {
        Self(ArkBn254::from(val))
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

impl PrimeField for Bn254 {
    fn modulus() -> BigUint {
        <ArkBn254 as ArkPrimeField>::MODULUS.into()
    }

    fn from_biguint(value: &BigUint) -> Self {
        Self(ArkBn254::from(value.clone()))
    }

    fn generator() -> BigUint {
        BigUint::from(5u32)
    }
}

impl PrimeFieldExt for Bn254 {
    fn to_biguint(&self) -> BigUint {
        self.0.into()
    }
}

impl PrimeFieldWords for Bn254 {
    fn to_words_le(&self) -> [u64; 4] {
        let limbs = self.0.into_bigint();
        let src = limbs.as_ref();
        let mut out = [0u64; 4];
        for (i, limb) in src.iter().enumerate().take(4) {
            out[i] = *limb;
        }
        out
    }
}
