use crate::fields::FieldElement;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};

pub trait MonolithField64: FieldElement {
    fn to_u64(&self) -> u64;
    fn modulus_u64() -> u64;
}

pub trait MonolithField32: FieldElement {
    fn to_u32(&self) -> u32;
    fn modulus_u32() -> u32;
}

#[derive(Clone, Debug)]
pub struct Monolith64Params<F: MonolithField64> {
    pub(crate) t: usize,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<Vec<F>>,
    pub(crate) mds: Vec<Vec<F>>,
    pub(crate) lookup: Vec<u16>,
}

#[derive(Clone, Debug)]
pub struct Monolith31Params<F: MonolithField32> {
    pub(crate) t: usize,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<Vec<F>>,
    pub(crate) mds: Vec<Vec<F>>,
    pub(crate) lookup1: Vec<u16>,
    pub(crate) lookup2: Vec<u16>,
}

impl<F: MonolithField64> Monolith64Params<F> {
    pub const R: usize = 6;
    pub const BARS: usize = 4;
    pub const INIT_SHAKE: &'static str = "Monolith";

    pub fn new(t: usize) -> Self {
        assert!(t == 8 || t == 12);
        let modulus = F::modulus_u64();
        let round_constants = instantiate_rc_64::<F>(t, modulus);
        let lookup = instantiate_lookup();
        let mds = mds_64::<F>(t);

        Monolith64Params {
            t,
            rounds: Self::R,
            round_constants,
            mds,
            lookup,
        }
    }

    pub fn get_t(&self) -> usize {
        self.t
    }

    pub fn get_rounds(&self) -> usize {
        self.rounds
    }
}

impl<F: MonolithField32> Monolith31Params<F> {
    pub const R: usize = 6;
    pub const BARS: usize = 8;
    pub const INIT_SHAKE: &'static str = "Monolith";

    pub fn new(t: usize) -> Self {
        assert!(t == 16 || t == 24);
        let modulus = F::modulus_u32();
        let round_constants = instantiate_rc_32::<F>(t, modulus);
        let lookup1 = instantiate_lookup();
        let lookup2 = instantiate_lookup2();
        let mds = mds_31::<F>(t);

        Monolith31Params {
            t,
            rounds: Self::R,
            round_constants,
            mds,
            lookup1,
            lookup2,
        }
    }

    pub fn get_t(&self) -> usize {
        self.t
    }

    pub fn get_rounds(&self) -> usize {
        self.rounds
    }
}

fn init_shake_64(t: usize, modulus: u64) -> impl XofReader {
    let mut shake = Shake128::default();
    shake.update(b"Monolith");
    shake.update(&[t as u8, 6u8]);
    shake.update(&u64::to_le_bytes(modulus));
    shake.update(&[8, 8, 8, 8, 8, 8, 8, 8]);
    shake.finalize_xof()
}

fn init_shake_32(t: usize, modulus: u32) -> impl XofReader {
    let mut shake = Shake128::default();
    shake.update(b"Monolith");
    shake.update(&[t as u8, 6u8]);
    shake.update(&u32::to_le_bytes(modulus));
    shake.update(&[8, 8, 8, 7]);
    shake.finalize_xof()
}

fn instantiate_rc_64<F: MonolithField64>(t: usize, modulus: u64) -> Vec<Vec<F>> {
    let mut reader = init_shake_64(t, modulus);
    (0..Monolith64Params::<F>::R - 1)
        .map(|_| {
            let mut row = Vec::with_capacity(t);
            for _ in 0..t {
                loop {
                    let mut buf = [0u8; 8];
                    reader.read(&mut buf);
                    let candidate = u64::from_le_bytes(buf);
                    if candidate < modulus {
                        row.push(F::from_u64(candidate));
                        break;
                    }
                }
            }
            row
        })
        .collect()
}

fn instantiate_rc_32<F: MonolithField32>(t: usize, modulus: u32) -> Vec<Vec<F>> {
    let mut reader = init_shake_32(t, modulus);
    (0..Monolith31Params::<F>::R - 1)
        .map(|_| {
            let mut row = Vec::with_capacity(t);
            for _ in 0..t {
                loop {
                    let mut buf = [0u8; 4];
                    reader.read(&mut buf);
                    let candidate = u32::from_le_bytes(buf);
                    if candidate < modulus {
                        row.push(F::from_u64(candidate as u64));
                        break;
                    }
                }
            }
            row
        })
        .collect()
}

fn bar0_8(limb: u8) -> u8 {
    let limbl1 = (limb >> 7) | (limb << 1);
    let limbl2 = (limb >> 6) | (limb << 2);
    let limbl3 = (limb >> 5) | (limb << 3);

    let tmp = limb ^ !limbl1 & limbl2 & limbl3;
    (tmp >> 7) | (tmp << 1)
}

fn bar1_7(limb: u8) -> u8 {
    let limbl1 = (limb >> 6) | (limb << 1);
    let limbl2 = (limb >> 5) | (limb << 2);

    let tmp = (limb ^ !limbl1 & limbl2) & 0x7f;
    ((tmp >> 6) | (tmp << 1)) & 0x7f
}

fn instantiate_lookup() -> Vec<u16> {
    (0..=u16::MAX)
        .map(|i| {
            let a = (i >> 8) as u8;
            let b = i as u8;
            ((bar0_8(a) as u16) << 8) | bar0_8(b) as u16
        })
        .collect()
}

fn instantiate_lookup2() -> Vec<u16> {
    (0..(1u32 << 15))
        .map(|i| {
            let a = (i >> 8) as u8;
            let b = i as u8;
            ((bar1_7(a) as u16) << 8) | bar0_8(b) as u16
        })
        .collect()
}

fn circulant_from_row<F: FieldElement>(row: &[u64], t: usize) -> Vec<Vec<F>> {
    assert!(row.len() >= t);
    let mut rot: Vec<F> = row.iter().map(|&v| F::from_u64(v)).collect();
    let mut full = Vec::with_capacity(row.len());
    full.push(rot.clone());
    for _ in 1..row.len() {
        rot.rotate_right(1);
        full.push(rot.clone());
    }

    full[..t].iter().map(|r| r[..t].to_vec()).collect()
}
// -- Circulant MDS matrix helpers --

fn mds_64<F: FieldElement>(t: usize) -> Vec<Vec<F>> {
    match t {
        8 => circulant_from_row::<F>(&[23, 8, 13, 10, 7, 6, 21, 8], t),
        12 => circulant_from_row::<F>(&[7, 23, 8, 26, 13, 10, 9, 7, 6, 22, 21, 8], t),
        _ => panic!("unsupported Monolith64 width"),
    }
}

fn mds_31<F: FieldElement>(t: usize) -> Vec<Vec<F>> {
    match t {
        16 => circulant_from_row::<F>(
            &[
                61402, 17845, 26798, 59689, 12021, 40901, 41351, 27521, 56951, 12034, 53865, 43244,
                7454, 33823, 28750, 1108,
            ],
            t,
        ),
        24 => circulant_from_row::<F>(
            &[
                87474966, 500304516, 1138910529, 1387408269, 937082352, 1410252806, 806711693,
                1520034124, 593719941, 1284124534, 1575767662, 927918294, 669885656, 1717383379,
                853820823, 1137173171, 1740948995, 2024301343, 1160738787, 60752863, 1950203872,
                1302354504, 1593997632, 136918578, 1358088042, 2071410473, 1467869360, 1941039814,
                1490713897, 1739211637, 230334003, 643163553,
            ],
            t,
        ),
        _ => panic!("unsupported Monolith31 width"),
    }
}
