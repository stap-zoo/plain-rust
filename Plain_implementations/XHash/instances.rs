use super::xhash::{InverseLayer, XHash, XHashParams, XHashProfile, XHASH_ROUNDS};
use crate::fields::goldilocks::Goldilocks;
use crate::fields::mersenne31::Mersenne31;
use crate::fields::FieldElement;
use lazy_static::lazy_static;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use std::sync::Arc;

const GOLDILOCKS_WIDTH: usize = 12;
const GOLDILOCKS_CAPACITY: usize = 4;
const GOLDILOCKS_ALPHA: u64 = 7;
const GOLDILOCKS_ALPHA_INV: u64 = 10_540_996_611_094_048_183;

const M31_WIDTH: usize = 24;
const M31_CAPACITY: usize = 8;
const M31_ALPHA: u64 = 5;
const M31_ALPHA_INV: u64 = 1_717_986_917;

/// First row of the 12x12 circulant RPO/XHash Goldilocks MDS matrix.
const GOLDILOCKS_MDS_FIRST_ROW: [u64; GOLDILOCKS_WIDTH] =
    [7, 23, 8, 26, 13, 10, 9, 7, 6, 22, 21, 8];

/// First row of the paper's 32x32 M31 circulant. The formal XHash M31
/// matrix is its top-left 24x24 submatrix, not a 24x24 circulant.
const M31_MDS_FIRST_ROW_32: [u32; 32] = [
    185_870_542,
    2_144_994_796,
    1_696_461_115,
    215_190_769,
    930_115_258,
    766_567_118,
    2_003_379_079,
    1_770_558_586,
    1_779_722_644,
    434_368_282,
    289_154_277,
    1_979_813_463,
    1_436_360_233,
    1_342_944_808,
    63_026_005,
    903_393_155,
    1_512_525_948,
    105_409_451,
    1_072_974_295,
    979_558_870,
    436_105_640,
    2_126_764_826,
    1_981_550_821,
    636_196_459,
    645_360_517,
    412_540_024,
    1_649_351_985,
    1_485_803_845,
    53_244_687,
    719_457_988,
    270_924_307,
    82_564_914,
];

lazy_static! {
    /// XHash8-Goldilocks: width 12 with eight active inverse S-boxes.
    pub static ref XHASH8_GOLDILOCKS_PARAMS: Arc<XHashParams<Goldilocks>> =
        Arc::new(goldilocks_params(InverseLayer::Partial));

    /// XHash12-Goldilocks: width 12 with a full inverse S-box layer.
    pub static ref XHASH12_GOLDILOCKS_PARAMS: Arc<XHashParams<Goldilocks>> =
        Arc::new(goldilocks_params(InverseLayer::Full));

    /// XHash16-M31: width 24 with sixteen active inverse S-boxes.
    pub static ref XHASH16_M31_PARAMS: Arc<XHashParams<Mersenne31>> =
        Arc::new(m31_params(InverseLayer::Partial));

    /// XHash24-M31: width 24 with a full inverse S-box layer.
    pub static ref XHASH24_M31_PARAMS: Arc<XHashParams<Mersenne31>> =
        Arc::new(m31_params(InverseLayer::Full));
}

pub fn xhash8_goldilocks() -> XHash<Goldilocks> {
    XHash::new(&XHASH8_GOLDILOCKS_PARAMS)
}

pub fn xhash12_goldilocks() -> XHash<Goldilocks> {
    XHash::new(&XHASH12_GOLDILOCKS_PARAMS)
}

pub fn xhash16_m31() -> XHash<Mersenne31> {
    XHash::new(&XHASH16_M31_PARAMS)
}

pub fn xhash24_m31() -> XHash<Mersenne31> {
    XHash::new(&XHASH24_M31_PARAMS)
}

fn goldilocks_params(inverse_layer: InverseLayer) -> XHashParams<Goldilocks> {
    let constants = generate_round_constants::<Goldilocks>(
        Goldilocks::MODULUS,
        GOLDILOCKS_WIDTH,
        GOLDILOCKS_CAPACITY,
    );
    XHashParams::new(
        GOLDILOCKS_WIDTH,
        GOLDILOCKS_CAPACITY,
        XHASH_ROUNDS,
        GOLDILOCKS_ALPHA,
        GOLDILOCKS_ALPHA_INV,
        XHashProfile::Goldilocks,
        inverse_layer,
        &constants,
        goldilocks_mds_matrix(),
    )
}

fn m31_params(inverse_layer: InverseLayer) -> XHashParams<Mersenne31> {
    let constants =
        generate_round_constants::<Mersenne31>(Mersenne31::MODULUS as u64, M31_WIDTH, M31_CAPACITY);
    XHashParams::new(
        M31_WIDTH,
        M31_CAPACITY,
        XHASH_ROUNDS,
        M31_ALPHA,
        M31_ALPHA_INV,
        XHashProfile::Mersenne31,
        inverse_layer,
        &constants,
        m31_mds_matrix(),
    )
}

/// Reproduce the official `XHash-Constants/round_constants.sage` generator
/// (commit 0045773de956d29d7246439d6c3dc6dce17c0586).
fn generate_round_constants<F: FieldElement>(
    modulus: u64,
    width: usize,
    capacity: usize,
) -> Vec<F> {
    let bit_length = (u64::BITS - modulus.leading_zeros()) as usize;
    let bytes_per_integer = bit_length.div_ceil(8) + 1;
    let constant_count = 3 * width * XHASH_ROUNDS;
    let seed = format!("XHash({modulus},{width},{capacity})");

    let mut hasher = Shake256::default();
    hasher.update(seed.as_bytes());
    let mut reader = hasher.finalize_xof();
    let mut bytes = vec![0u8; bytes_per_integer * constant_count];
    reader.read(&mut bytes);

    bytes
        .chunks_exact(bytes_per_integer)
        .map(|chunk| {
            let integer = chunk.iter().enumerate().fold(0u128, |acc, (index, &byte)| {
                acc + ((byte as u128) << (8 * index))
            });
            F::from_u64((integer % modulus as u128) as u64)
        })
        .collect()
}

fn goldilocks_mds_matrix() -> Vec<Vec<Goldilocks>> {
    (0..GOLDILOCKS_WIDTH)
        .map(|row| {
            (0..GOLDILOCKS_WIDTH)
                .map(|column| {
                    let index = (column + GOLDILOCKS_WIDTH - row) % GOLDILOCKS_WIDTH;
                    Goldilocks::from_u64(GOLDILOCKS_MDS_FIRST_ROW[index])
                })
                .collect()
        })
        .collect()
}

fn m31_mds_matrix() -> Vec<Vec<Mersenne31>> {
    (0..M31_WIDTH)
        .map(|row| {
            (0..M31_WIDTH)
                .map(|column| {
                    let index = (column + 32 - row) % 32;
                    Mersenne31::from_u32(M31_MDS_FIRST_ROW_32[index])
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{
        M31_MDS_FIRST_ROW_32, XHASH12_GOLDILOCKS_PARAMS, XHASH16_M31_PARAMS, XHASH24_M31_PARAMS,
        XHASH8_GOLDILOCKS_PARAMS,
    };

    #[test]
    fn official_constant_prefixes_and_lengths_match() {
        let gold = &XHASH8_GOLDILOCKS_PARAMS.round_constants;
        assert_eq!(gold.len(), 3 * 12 * 6);
        assert_eq!(gold[0].to_u64(), 2_342_829_598_172_340_146);
        assert_eq!(gold[1].to_u64(), 3_585_263_893_330_295_736);
        assert_eq!(gold[2].to_u64(), 9_840_230_240_082_666_977);
        assert_eq!(gold[gold.len() - 3].to_u64(), 14_626_303_366_643_264_814);
        assert_eq!(gold[gold.len() - 2].to_u64(), 6_015_835_057_058_442_672);
        assert_eq!(gold[gold.len() - 1].to_u64(), 3_697_516_607_449_911_960);

        let m31 = &XHASH16_M31_PARAMS.round_constants;
        assert_eq!(m31.len(), 3 * 24 * 6);
        assert_eq!(m31[0].to_u32(), 175_084_324);
        assert_eq!(m31[1].to_u32(), 307_267_372);
        assert_eq!(m31[2].to_u32(), 926_126_032);
        assert_eq!(m31[m31.len() - 3].to_u32(), 1_156_311_902);
        assert_eq!(m31[m31.len() - 2].to_u32(), 225_372_763);
        assert_eq!(m31[m31.len() - 1].to_u32(), 1_906_892_131);

        assert_eq!(
            XHASH8_GOLDILOCKS_PARAMS.round_constants,
            XHASH12_GOLDILOCKS_PARAMS.round_constants
        );
        assert_eq!(
            XHASH16_M31_PARAMS.round_constants,
            XHASH24_M31_PARAMS.round_constants
        );
    }

    #[test]
    fn m31_matrix_is_the_full_32_circulant_truncation() {
        let matrix = XHASH24_M31_PARAMS.mds_matrix();
        assert_eq!(matrix[0][0].to_u32(), M31_MDS_FIRST_ROW_32[0]);
        assert_eq!(matrix[0][23].to_u32(), M31_MDS_FIRST_ROW_32[23]);
        assert_eq!(matrix[1][0].to_u32(), M31_MDS_FIRST_ROW_32[31]);
        assert_eq!(matrix[23][0].to_u32(), M31_MDS_FIRST_ROW_32[9]);
    }
}
