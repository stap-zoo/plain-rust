use super::monolith_params::{Monolith31Params, MonolithField32};
use crate::fields::koalabear::KoalaBear;
use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = KoalaBear;

impl MonolithField32 for KoalaBear {
    fn to_u32(&self) -> u32 {
        KoalaBear::to_u32(self)
    }

    fn modulus_u32() -> u32 {
        KoalaBear::MODULUS
    }
}

lazy_static! {
    pub static ref MONOLITH_KOALABEAR_16_PARAMS: Arc<Monolith31Params<Scalar>> =
        Arc::new(Monolith31Params::new(16));
    pub static ref MONOLITH_KOALABEAR_24_PARAMS: Arc<Monolith31Params<Scalar>> =
        Arc::new(Monolith31Params::new(24));
}
