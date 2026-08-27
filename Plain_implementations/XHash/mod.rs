pub mod instances;
pub mod xhash;

pub use instances::{
    xhash12_goldilocks, xhash16_m31, xhash24_m31, xhash8_goldilocks, XHASH12_GOLDILOCKS_PARAMS,
    XHASH16_M31_PARAMS, XHASH24_M31_PARAMS, XHASH8_GOLDILOCKS_PARAMS,
};
pub use xhash::{InverseLayer, XHash, XHashParams, XHashProfile, XHASH_ROUNDS};
