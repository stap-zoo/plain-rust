use super::monolith_params::{
    Monolith31Params, Monolith64Params, MonolithField32, MonolithField64,
};
use crate::fields::babybear::BabyBear;
use crate::fields::goldilocks::Goldilocks;
use crate::fields::koalabear::KoalaBear;
use crate::fields::mersenne31::Mersenne31;
use lazy_static::lazy_static;
use std::sync::Arc;

// ============================================================
// MonolithField implementations
// ============================================================

impl MonolithField64 for Goldilocks {
    fn to_u64(&self) -> u64 {
        Goldilocks::to_u64(self)
    }

    fn modulus_u64() -> u64 {
        Goldilocks::MODULUS
    }
}

impl MonolithField32 for Mersenne31 {
    fn to_u32(&self) -> u32 {
        Mersenne31::to_u32(self)
    }

    fn modulus_u32() -> u32 {
        Mersenne31::MODULUS
    }
}

impl MonolithField32 for BabyBear {
    fn to_u32(&self) -> u32 {
        BabyBear::to_u32(self)
    }

    fn modulus_u32() -> u32 {
        BabyBear::MODULUS
    }
}

impl MonolithField32 for KoalaBear {
    fn to_u32(&self) -> u32 {
        KoalaBear::to_u32(self)
    }

    fn modulus_u32() -> u32 {
        KoalaBear::MODULUS
    }
}

lazy_static! {
    pub static ref MONOLITH_GOLDILOCKS_8_PARAMS: Arc<Monolith64Params<Goldilocks>> =
        Arc::new(Monolith64Params::new(8));
    pub static ref MONOLITH_GOLDILOCKS_12_PARAMS: Arc<Monolith64Params<Goldilocks>> =
        Arc::new(Monolith64Params::new(12));
    pub static ref MONOLITH_MERSENNE31_16_PARAMS: Arc<Monolith31Params<Mersenne31>> =
        Arc::new(Monolith31Params::new(16));
    pub static ref MONOLITH_MERSENNE31_24_PARAMS: Arc<Monolith31Params<Mersenne31>> =
        Arc::new(Monolith31Params::new(24));
    pub static ref MONOLITH_BABYBEAR_16_PARAMS: Arc<Monolith31Params<BabyBear>> =
        Arc::new(Monolith31Params::new(16));
    pub static ref MONOLITH_BABYBEAR_24_PARAMS: Arc<Monolith31Params<BabyBear>> =
        Arc::new(Monolith31Params::new(24));
    pub static ref MONOLITH_KOALABEAR_16_PARAMS: Arc<Monolith31Params<KoalaBear>> =
        Arc::new(Monolith31Params::new(16));
    pub static ref MONOLITH_KOALABEAR_24_PARAMS: Arc<Monolith31Params<KoalaBear>> =
        Arc::new(Monolith31Params::new(24));
}
