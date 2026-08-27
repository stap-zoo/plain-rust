use super::{FieldElement, PrimeField, PrimeFieldExt, PrimeFieldWords};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use p3_field::{Field as P3Field, PrimeCharacteristicRing, PrimeField64};
use p3_goldilocks::Goldilocks as P3Goldilocks;

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct Goldilocks(pub(crate) P3Goldilocks);

impl Goldilocks {
    pub const MODULUS: u64 = <P3Goldilocks as PrimeField64>::ORDER_U64;

    pub fn from_hex(s: &str) -> Option<Self> {
        let s = s.strip_prefix("0x").unwrap_or(s);
        u64::from_str_radix(s, 16).ok().map(Self::from_u64)
    }

    pub fn to_u64(&self) -> u64 {
        <P3Goldilocks as PrimeField64>::as_canonical_u64(&self.0)
    }
}

impl FieldElement for Goldilocks {
    fn zero() -> Self {
        Self(P3Goldilocks::from_u64(0))
    }

    fn one() -> Self {
        Self(P3Goldilocks::from_u64(1))
    }

    fn from_u64(val: u64) -> Self {
        Self(P3Goldilocks::from_u64(val))
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

impl PrimeField for Goldilocks {
    fn modulus() -> BigUint {
        <P3Goldilocks as P3Field>::order()
    }

    fn from_biguint(value: &BigUint) -> Self {
        let modulus = Self::modulus();
        let reduced = value % &modulus;
        let v = reduced.to_u64().expect("Goldilocks fits into u64");
        Self::from_u64(v)
    }

    fn generator() -> BigUint {
        BigUint::from(7u32)
    }
}

impl PrimeFieldExt for Goldilocks {
    fn to_biguint(&self) -> BigUint {
        BigUint::from(self.to_u64())
    }
}

impl PrimeFieldWords for Goldilocks {
    fn to_words_le(&self) -> [u64; 4] {
        [self.to_u64(), 0, 0, 0]
    }
}
