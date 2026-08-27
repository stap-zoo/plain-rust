use super::monolith_params::{Monolith64Params, MonolithField64};
use crate::fields::goldilocks::Goldilocks;
use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = Goldilocks;

impl MonolithField64 for Goldilocks {
    fn to_u64(&self) -> u64 {
        Goldilocks::to_u64(self)
    }

    fn modulus_u64() -> u64 {
        Goldilocks::MODULUS
    }
}

lazy_static! {
    pub static ref MONOLITH_GOLDILOCKS_8_PARAMS: Arc<Monolith64Params<Scalar>> =
        Arc::new(Monolith64Params::new(8));
    pub static ref MONOLITH_GOLDILOCKS_12_PARAMS: Arc<Monolith64Params<Scalar>> =
        Arc::new(Monolith64Params::new(12));
}
