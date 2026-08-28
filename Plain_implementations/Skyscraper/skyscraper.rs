use crate::fields::{PrimeField, PrimeFieldMontgomery};
use std::sync::Arc;

const TOTAL_ROUNDS: usize = 18;
const BAR_ROUNDS: [usize; 4] = [6, 7, 10, 11];

#[derive(Clone, Debug)]
pub struct ExtElem<F: PrimeField, const N: usize> {
    coeffs: [F; N],
}

impl<F: PrimeField, const N: usize> ExtElem<F, N> {
    fn from_coeffs(coeffs: [F; N]) -> Self {
        ExtElem { coeffs }
    }

    fn coeffs(&self) -> &[F; N] {
        &self.coeffs
    }

    fn add_assign(&mut self, other: &Self) {
        for (a, b) in self.coeffs.iter_mut().zip(other.coeffs.iter()) {
            a.add_assign(b);
        }
    }
}

/// The Square and Bar sboxes below operate on coefficients kept in *raw* Montgomery-domain
/// form throughout the permutation (see `PrimeFieldMontgomery`): every field multiply on such
/// values carries its own free `* R^-1` reduction, which is exactly Skyscraper's specified
/// `(a+bX)^2 mod (X^2+beta)` sbox for the Square rounds -- no separate correction multiply
/// needed -- and it lets the Bar rounds read/write a coefficient's canonical integer directly,
/// with no Montgomery conversion at all. `beta` and the round constants are prepared to match
/// this convention in `SkyscraperParams::new` and `Skyscraper::permutation`.
impl<F: PrimeFieldMontgomery, const N: usize> ExtElem<F, N> {
    fn square_in_place(&mut self, beta: &F) {
        match N {
            2 => self.square_n2(beta),
            3 => self.square_n3(beta),
            _ => panic!("unsupported extension degree"),
        }
    }

    fn square_n2(&mut self, beta: &F) {
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

        self.coeffs[0] = out0;
        self.coeffs[1] = out1;
    }

    fn square_n3(&mut self, beta: &F) {
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

        self.coeffs[0] = out0;
        self.coeffs[1] = out1;
        self.coeffs[2] = out2;
    }

    fn bar_in_place(&mut self) {
        let mut lows = [0u128; N];
        let mut highs = [0u128; N];

        for idx in 0..N {
            let (low, high) = words_to_halves(self.coeffs[idx].raw_words());
            lows[idx] = low;
            highs[idx] = high;
        }

        for idx in 0..N {
            let prev = (idx + N - 1) % N;
            let low = bar_u128(highs[prev]);
            let high = bar_u128(lows[idx]);
            self.coeffs[idx] = F::from_reduced_raw_words(halves_to_words(low, high));
        }
    }
}

#[derive(Clone, Debug)]
pub struct SkyscraperParams<F: PrimeField, const N: usize> {
    pub(crate) beta: F,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<ExtElem<F, N>>,
}

impl<F: PrimeFieldMontgomery, const N: usize> SkyscraperParams<F, N> {
    pub fn new(beta: u64, round_constants: &[Vec<F>]) -> Self {
        assert!(N == 2 || N == 3);
        assert_eq!(round_constants.len(), TOTAL_ROUNDS - 2);
        for rc in round_constants {
            assert_eq!(rc.len(), N);
        }

        // The old static tables were laid out high-coordinate first. The reference
        // state uses polynomial coefficients in ascending degree order. Round constants
        // are converted into the same raw Montgomery-domain form the permutation state
        // uses, so they can be added directly into it every round.
        let round_constants = round_constants
            .iter()
            .map(|row| -> Vec<F> {
                row.iter()
                    .rev()
                    .map(|c| c.clone().into_montgomery_raw())
                    .collect()
            })
            .map(|coeffs| {
                ExtElem::from_coeffs(
                    coeffs
                        .try_into()
                        .unwrap_or_else(|_| panic!("expected {N} round-constant coefficients")),
                )
            })
            .collect();

        SkyscraperParams {
            beta: F::from_u64(beta),
            rounds: TOTAL_ROUNDS,
            round_constants,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Skyscraper<F: PrimeFieldMontgomery, const N: usize> {
    pub(crate) params: Arc<SkyscraperParams<F, N>>,
}

impl<F: PrimeFieldMontgomery, const N: usize> Skyscraper<F, N> {
    pub fn new(params: &Arc<SkyscraperParams<F, N>>) -> Self {
        Skyscraper {
            params: Arc::clone(params),
        }
    }

    pub fn get_n(&self) -> usize {
        N
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        assert_eq!(input.len(), 2 * N);

        let mut left = ExtElem::from_coeffs(std::array::from_fn(|i| {
            input[i].clone().into_montgomery_raw()
        }));
        let mut right = ExtElem::from_coeffs(std::array::from_fn(|i| {
            input[N + i].clone().into_montgomery_raw()
        }));

        for round in 0..self.params.rounds {
            let prev_left = left.clone();

            if is_bar_round(round) {
                left.bar_in_place();
            } else {
                left.square_in_place(&self.params.beta);
            }

            if (1..(self.params.rounds - 1)).contains(&round) {
                right.add_assign(&self.params.round_constants[round - 1]);
            }

            left.add_assign(&right);
            right = prev_left;
        }

        let mut out = Vec::with_capacity(2 * N);
        out.extend(left.coeffs().iter().cloned().map(F::from_montgomery_raw));
        out.extend(right.coeffs().iter().cloned().map(F::from_montgomery_raw));
        out
    }
}

#[inline(always)]
fn is_bar_round(round: usize) -> bool {
    BAR_ROUNDS.contains(&round)
}

fn words_to_halves(words: [u64; 4]) -> (u128, u128) {
    let low = u128::from(words[0]) | (u128::from(words[1]) << 64);
    let high = u128::from(words[2]) | (u128::from(words[3]) << 64);
    (low, high)
}

fn halves_to_words(low: u128, high: u128) -> [u64; 4] {
    [
        low as u64,
        (low >> 64) as u64,
        high as u64,
        (high >> 64) as u64,
    ]
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
