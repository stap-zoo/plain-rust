use super::monolith_params::{Monolith31Params, MonolithField32};
use crate::fields::mersenne31::Mersenne31;
use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = Mersenne31;

impl MonolithField32 for Mersenne31 {
    fn to_u32(&self) -> u32 {
        Mersenne31::to_u32(self)
    }

    fn modulus_u32() -> u32 {
        Mersenne31::MODULUS
    }
}

lazy_static! {
    pub static ref MONOLITH_MERSENNE31_16_PARAMS: Arc<Monolith31Params<Scalar>> =
        Arc::new(Monolith31Params::new(16));
    pub static ref MONOLITH_MERSENNE31_24_PARAMS: Arc<Monolith31Params<Scalar>> =
        Arc::new(Monolith31Params::new(24));
}
