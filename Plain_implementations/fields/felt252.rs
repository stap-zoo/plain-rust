use super::montgomery_4::{
    add_mod, from_hex_to_limbs, monty_mul, reduce_raw, sub_mod, to_monty, MontyParams,
};
use super::{biguint_from_limbs_le, FieldElement, PrimeField, PrimeFieldExt, PrimeFieldWords};
use num_bigint::BigUint;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Felt252 {
    pub(crate) value: [u64; 4],
}

struct Felt252Params;

impl MontyParams for Felt252Params {
    const MODULUS: [u64; 4] = [0x1, 0x0, 0x0, 0x0800000000000011];
    const MU: u64 = 0x1;
    const R2: [u64; 4] = [
        0xfffffd737e000401,
        0x00000001330fffff,
        0xffffffffff6f8000,
        0x07ffd4ab5e008810,
    ];
}

impl Felt252 {
    pub fn from_hex(s: &str) -> Option<Self> {
        let raw = from_hex_to_limbs(s)?;
        Some(Self::from_raw(raw))
    }

    #[inline]
    fn from_raw(raw: [u64; 4]) -> Self {
        Self {
            value: to_monty::<Felt252Params>(reduce_raw::<Felt252Params>(raw)),
        }
    }
}

impl FieldElement for Felt252 {
    fn zero() -> Self {
        Self { value: [0; 4] }
    }

    fn one() -> Self {
        Self::from_u64(1)
    }

    fn from_u64(val: u64) -> Self {
        Self::from_raw([val, 0, 0, 0])
    }

    fn add_assign(&mut self, other: &Self) {
        self.value = add_mod::<Felt252Params>(self.value, other.value);
    }

    fn sub_assign(&mut self, other: &Self) {
        self.value = sub_mod::<Felt252Params>(self.value, other.value);
    }

    fn mul_assign(&mut self, other: &Self) {
        self.value = monty_mul::<Felt252Params>(self.value, other.value);
    }
}

impl PrimeField for Felt252 {
    fn modulus() -> BigUint {
        BigUint::parse_bytes(
            b"800000000000011000000000000000000000000000000000000000000000001",
            16,
        )
        .expect("valid felt252 modulus")
    }

    fn from_biguint(value: &BigUint) -> Self {
        let modulus = Self::modulus();
        let reduced = value % &modulus;
        let hex = reduced.to_str_radix(16);
        let hex = if hex.is_empty() { "0".to_string() } else { hex };
        let prefixed = format!("0x{hex}");
        Self::from_hex(&prefixed).expect("valid felt252 element")
    }

    fn generator() -> BigUint {
        BigUint::from(3u32)
    }
}

impl PrimeFieldExt for Felt252 {
    fn to_biguint(&self) -> BigUint {
        let normal = monty_mul::<Felt252Params>(self.value, [1, 0, 0, 0]);
        biguint_from_limbs_le(&normal)
    }
}

impl PrimeFieldWords for Felt252 {
    fn to_words_le(&self) -> [u64; 4] {
        monty_mul::<Felt252Params>(self.value, [1, 0, 0, 0])
    }
}
