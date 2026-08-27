use super::monolith_params::{Monolith31Params, MonolithField32};
use crate::fields::babybear::BabyBear;
use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = BabyBear;

impl MonolithField32 for BabyBear {
    fn to_u32(&self) -> u32 {
        BabyBear::to_u32(self)
    }

    fn modulus_u32() -> u32 {
        BabyBear::MODULUS
    }
}

lazy_static! {
    pub static ref MONOLITH_BABYBEAR_16_PARAMS: Arc<Monolith31Params<Scalar>> =
        Arc::new(Monolith31Params::new(16));
    pub static ref MONOLITH_BABYBEAR_24_PARAMS: Arc<Monolith31Params<Scalar>> =
        Arc::new(Monolith31Params::new(24));
}
