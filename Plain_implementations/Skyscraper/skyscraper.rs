use crate::fields::{PrimeField, PrimeFieldWords};
use num_bigint::BigUint;
use num_traits::One;
use std::sync::Arc;

const TOTAL_ROUNDS: usize = 18;
const BAR_ROUNDS: [usize; 4] = [6, 7, 10, 11];

#[derive(Clone, Debug)]
pub struct ExtElem<F: PrimeField> {
    coeffs: Vec<F>,
}

impl<F: PrimeField> ExtElem<F> {
    fn from_coeffs(coeffs: Vec<F>) -> Self {
        ExtElem { coeffs }
    }

    fn add_assign(&mut self, other: &Self) {
        for (a, b) in self.coeffs.iter_mut().zip(other.coeffs.iter()) {
            a.add_assign(b);
        }
    }

    fn square_in_place(&mut self, beta: &F, mont_r_inv: &F) {
        match self.coeffs.len() {
            2 => self.square_n2(beta, mont_r_inv),
            3 => self.square_n3(beta, mont_r_inv),
            _ => panic!("unsupported extension degree"),
        }
    }

    fn square_n2(&mut self, beta: &F, mont_r_inv: &F) {
        let a = self.coeffs[0].clone();
        let b = self.coeffs[1].clone();

        // (a + bX)^2 mod (X^2 + beta).
        let mut out0 = a.clone();
        out0.square();
        let mut beta_b2 = b.clone();
        beta_b2.square();
        beta_b2.mul_assign(beta);
        out0.sub_assign(&beta_b2);

        let mut out1 = a;
        out1.mul_assign(&b);
        out1.double();

        out0.mul_assign(mont_r_inv);
        out1.mul_assign(mont_r_inv);

        self.coeffs[0] = out0;
        self.coeffs[1] = out1;
    }

    fn square_n3(&mut self, beta: &F, mont_r_inv: &F) {
        let a = self.coeffs[0].clone();
        let b = self.coeffs[1].clone();
        let c = self.coeffs[2].clone();

        // (a + bX + cX^2)^2 mod (X^3 + beta).
        let mut out0 = a.clone();
        out0.square();
        let mut two_beta_bc = b.clone();
        two_beta_bc.mul_assign(&c);
        two_beta_bc.double();
        two_beta_bc.mul_assign(beta);
        out0.sub_assign(&two_beta_bc);

        let mut out1 = a.clone();
        out1.mul_assign(&b);
        out1.double();
        let mut beta_c2 = c.clone();
        beta_c2.square();
        beta_c2.mul_assign(beta);
        out1.sub_assign(&beta_c2);

        let mut out2 = b;
        out2.square();
        let mut two_ac = a;
        two_ac.mul_assign(&c);
        two_ac.double();
        out2.add_assign(&two_ac);

        out0.mul_assign(mont_r_inv);
        out1.mul_assign(mont_r_inv);
        out2.mul_assign(mont_r_inv);

        self.coeffs[0] = out0;
        self.coeffs[1] = out1;
        self.coeffs[2] = out2;
    }
}

impl<F: PrimeFieldWords> ExtElem<F> {
    fn bar_in_place(&mut self) {
        let n = self.coeffs.len();
        let mut lows = vec![0u128; n];
        let mut highs = vec![0u128; n];

        for (idx, coeff) in self.coeffs.iter().enumerate() {
            let (low, high) = split_halves_u128(coeff);
            lows[idx] = low;
            highs[idx] = high;
        }

        for idx in 0..n {
            let prev = (idx + n - 1) % n;
            let low = bar_u128(highs[prev]);
            let high = bar_u128(lows[idx]);
            self.coeffs[idx] = field_from_halves(low, high);
        }
    }
}

#[derive(Clone, Debug)]
pub struct SkyscraperParams<F: PrimeField> {
    pub(crate) n: usize,
    pub(crate) beta_f: F,
    pub(crate) mont_r_inv: F,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<ExtElem<F>>,
}

impl<F: PrimeField> SkyscraperParams<F> {
    pub fn new(n: usize, beta: u64, round_constants: &[Vec<F>]) -> Self {
        assert!(n == 2 || n == 3);
        assert_eq!(round_constants.len(), TOTAL_ROUNDS - 2);
        for rc in round_constants {
            assert_eq!(rc.len(), n);
        }

        // The old static tables were laid out high-coordinate first. The reference
        // state uses polynomial coefficients in ascending degree order.
        let round_constants = round_constants
            .iter()
            .map(|row| row.iter().rev().cloned().collect())
            .map(ExtElem::from_coeffs)
            .collect();

        let modulus = F::modulus();
        let machine_bits = (modulus.bits() as usize).div_ceil(64) * 64;
        let mont_r = F::from_biguint(&(BigUint::one() << machine_bits));
        let mont_r_inv = mont_r.pow_words_le(&(modulus - BigUint::from(2u64)).to_u64_digits());

        SkyscraperParams {
            n,
            beta_f: F::from_u64(beta),
            mont_r_inv,
            rounds: TOTAL_ROUNDS,
            round_constants,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Skyscraper<F: PrimeFieldWords> {
    pub(crate) params: Arc<SkyscraperParams<F>>,
}

impl<F: PrimeFieldWords> Skyscraper<F> {
    pub fn new(params: &Arc<SkyscraperParams<F>>) -> Self {
        Skyscraper {
            params: Arc::clone(params),
        }
    }

    pub fn get_n(&self) -> usize {
        self.params.n
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let n = self.params.n;
        assert_eq!(input.len(), 2 * n);

        let mut left = ExtElem::from_coeffs(input[..n].to_vec());
        let mut right = ExtElem::from_coeffs(input[n..].to_vec());

        for round in 0..self.params.rounds {
            let prev_left = left.clone();

            if is_bar_round(round) {
                left.bar_in_place();
            } else {
                left.square_in_place(&self.params.beta_f, &self.params.mont_r_inv);
            }

            if (1..(self.params.rounds - 1)).contains(&round) {
                right.add_assign(&self.params.round_constants[round - 1]);
            }

            left.add_assign(&right);
            right = prev_left;
        }

        let mut out = Vec::with_capacity(2 * n);
        out.extend_from_slice(&left.coeffs);
        out.extend_from_slice(&right.coeffs);
        out
    }
}

#[inline(always)]
fn is_bar_round(round: usize) -> bool {
    BAR_ROUNDS.contains(&round)
}

fn split_halves_u128<F: PrimeFieldWords>(value: &F) -> (u128, u128) {
    let words = value.to_words_le();
    let low = u128::from(words[0]) | (u128::from(words[1]) << 64);
    let high = u128::from(words[2]) | (u128::from(words[3]) << 64);
    (low, high)
}

fn field_from_halves<F: PrimeFieldWords>(low: u128, high: u128) -> F {
    F::from_words_le([
        low as u64,
        (low >> 64) as u64,
        high as u64,
        (high >> 64) as u64,
    ])
}

fn bar_u128(value: u128) -> u128 {
    const MASK_80: u128 = 0x80808080808080808080808080808080;
    const MASK_7F: u128 = 0x7f7f7f7f7f7f7f7f7f7f7f7f7f7f7f7f;
    const MASK_C0: u128 = 0xc0c0c0c0c0c0c0c0c0c0c0c0c0c0c0c0;
    const MASK_3F: u128 = 0x3f3f3f3f3f3f3f3f3f3f3f3f3f3f3f3f;
    const MASK_E0: u128 = 0xe0e0e0e0e0e0e0e0e0e0e0e0e0e0e0e0;
    const MASK_1F: u128 = 0x1f1f1f1f1f1f1f1f1f1f1f1f1f1f1f1f;

    let t1 = ((value & MASK_80) >> 7) | ((value & MASK_7F) << 1);
    let t2 = ((value & MASK_C0) >> 6) | ((value & MASK_3F) << 2);
    let t3 = ((value & MASK_E0) >> 5) | ((value & MASK_1F) << 3);
    let tmp = (!t1 & t2 & t3) ^ value;
    ((tmp & MASK_80) >> 7) | ((tmp & MASK_7F) << 1)
}
