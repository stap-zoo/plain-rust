use super::{FieldElement, PrimeField, PrimeFieldExt, PrimeFieldMontgomery, PrimeFieldWords};
use ark_bn254::Fr as ArkBn254;
use ark_ff::{BigInt, BigInteger, PrimeField as ArkPrimeField};
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

    fn from_words_le(words: [u64; 4]) -> Self {
        let mut bytes = [0u8; 32];
        for (chunk, word) in bytes.chunks_exact_mut(8).zip(words.iter()) {
            chunk.copy_from_slice(&word.to_le_bytes());
        }
        Self(ArkBn254::from_le_bytes_mod_order(&bytes))
    }
}

impl PrimeFieldMontgomery for Bn254 {
    fn into_montgomery_raw(self) -> Self {
        Self(ArkBn254::new_unchecked(self.0.into_bigint()))
    }

    fn from_montgomery_raw(self) -> Self {
        Self(ArkBn254::new(self.0 .0))
    }

    fn raw_words(&self) -> [u64; 4] {
        (self.0 .0).0
    }

    fn from_reduced_raw_words(words: [u64; 4]) -> Self {
        let mut value = BigInt(words);
        let modulus = <ArkBn254 as ArkPrimeField>::MODULUS;
        while value >= modulus {
            value.sub_with_borrow(&modulus);
        }
        Self(ArkBn254::new_unchecked(value))
    }
}

#[cfg(test)]
mod montgomery_raw_tests {
    use super::*;
    use num_traits::One;

    #[test]
    fn raw_square_applies_exactly_one_r_inverse() {
        let x = Bn254::from_u64(5);
        let mut x_raw = x.into_montgomery_raw();
        x_raw.square();
        let actual = x_raw.from_montgomery_raw();

        let modulus = Bn254::modulus();
        let machine_bits = (modulus.bits() as usize).div_ceil(64) * 64;
        let mont_r = Bn254::from_biguint(&(BigUint::one() << machine_bits));
        let mont_r_inv = mont_r.pow_words_le(&(modulus - BigUint::from(2u64)).to_u64_digits());

        let mut expected = Bn254::from_u64(25);
        expected.mul_assign(&mont_r_inv);

        assert_eq!(actual, expected, "x_raw.square().from_montgomery_raw() should equal x^2 * R^-1");
        assert_ne!(actual, Bn254::from_u64(25), "sanity: the raw trick must NOT equal the plain square");
    }

    #[test]
    fn raw_then_normal_round_trip_is_identity() {
        let x = Bn254::from_u64(1234567);
        let round_tripped = x.into_montgomery_raw().from_montgomery_raw();
        assert_eq!(x, round_tripped);
    }

    #[test]
    fn raw_words_reduction_matches_from_words_le() {
        // A value already < modulus should round-trip identically through both paths.
        let words = [0x1122334455667788u64, 0x99aabbccddeeff00, 0, 0];
        let via_raw = Bn254::from_reduced_raw_words(words).from_montgomery_raw();
        let via_normal = Bn254::from_words_le(words);
        assert_eq!(via_raw, via_normal);
    }
}
