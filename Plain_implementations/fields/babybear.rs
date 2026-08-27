use super::{FieldElement, PrimeField, PrimeFieldExt, PrimeFieldWords};
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use p3_baby_bear::BabyBear as P3BabyBear;
use p3_field::{Field as P3Field, PrimeCharacteristicRing, PrimeField32};

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct BabyBear(pub(crate) P3BabyBear);

impl BabyBear {
    pub const MODULUS: u32 = <P3BabyBear as PrimeField32>::ORDER_U32;

    pub fn from_u32(val: u32) -> Self {
        Self(P3BabyBear::from_u32(val))
    }

    pub fn to_u32(&self) -> u32 {
        <P3BabyBear as PrimeField32>::as_canonical_u32(&self.0)
    }
}

impl FieldElement for BabyBear {
    fn zero() -> Self {
        Self(P3BabyBear::from_u64(0))
    }

    fn one() -> Self {
        Self(P3BabyBear::from_u64(1))
    }

    fn from_u64(val: u64) -> Self {
        Self(P3BabyBear::from_u64(val))
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

impl PrimeField for BabyBear {
    fn modulus() -> BigUint {
        <P3BabyBear as P3Field>::order()
    }

    fn from_biguint(value: &BigUint) -> Self {
        let modulus = Self::modulus();
        let reduced = value % &modulus;
        let v = reduced.to_u64().expect("BabyBear fits into u64");
        Self::from_u32(v as u32)
    }

    fn generator() -> BigUint {
        BigUint::from(31u32)
    }
}

impl PrimeFieldExt for BabyBear {
    fn to_biguint(&self) -> BigUint {
        BigUint::from(self.to_u32() as u64)
    }
}

impl PrimeFieldWords for BabyBear {
    fn to_words_le(&self) -> [u64; 4] {
        [self.to_u32() as u64, 0, 0, 0]
    }
}
