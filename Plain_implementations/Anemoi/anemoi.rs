use crate::fields::FieldElement;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct AnemoiParams<F: FieldElement> {
    pub(crate) n_cols: usize,
    pub(crate) width: usize,
    pub(crate) rounds: usize,
    pub(crate) alpha: u64,
    pub(crate) alpha_inv: [u64; 4],
    pub(crate) beta: F,
    pub(crate) delta: F,
    pub(crate) mds: Vec<Vec<F>>,
    pub(crate) round_constants_c: Vec<Vec<F>>,
    pub(crate) round_constants_d: Vec<Vec<F>>,
}

#[derive(Clone, Debug)]
pub struct Anemoi<F: FieldElement> {
    pub(crate) params: Arc<AnemoiParams<F>>,
}

impl<F: FieldElement> Anemoi<F> {
    pub fn new(params: &Arc<AnemoiParams<F>>) -> Self {
        Anemoi {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.width
    }

    pub fn get_alpha(&self) -> u64 {
        self.params.alpha
    }

    pub fn get_rounds(&self) -> usize {
        self.params.rounds
    }

    pub fn permutation(&self, _input: &[F]) -> Vec<F> {
        let width = self.params.width;
        let n_cols = self.params.n_cols;
        assert_eq!(width, 2 * n_cols);
        assert_eq!(_input.len(), width);

        let mut state = _input.to_vec();
        for r in 0..self.params.rounds {
            self.add_round_constants(&mut state, r);
            self.linear_layer(&mut state);
            self.sbox_layer(&mut state);
        }
        self.linear_layer(&mut state);
        state
    }

    fn add_round_constants(&self, state: &mut [F], round: usize) {
        let n_cols = self.params.n_cols;
        for i in 0..n_cols {
            let mut x = state[i].clone();
            x.add_assign(&self.params.round_constants_c[round][i]);
            state[i] = x;
            let mut y = state[n_cols + i].clone();
            y.add_assign(&self.params.round_constants_d[round][i]);
            state[n_cols + i] = y;
        }
    }

    fn mds_mul(&self, input: &[F]) -> Vec<F> {
        let n_cols = self.params.n_cols;
        let mut output = vec![F::zero(); n_cols];
        for r in 0..n_cols {
            let mut acc = F::zero();
            for c in 0..n_cols {
                let mut tmp = self.params.mds[r][c].clone();
                tmp.mul_assign(&input[c]);
                acc.add_assign(&tmp);
            }
            output[r] = acc;
        }
        output
    }

    fn linear_layer(&self, state: &mut [F]) {
        let n_cols = self.params.n_cols;

        self.apply_mds_only(state);

        for i in 0..n_cols {
            let mut new_y = state[n_cols + i].clone();
            new_y.add_assign(&state[i]);
            state[n_cols + i] = new_y;

            let mut new_x = state[i].clone();
            new_x.add_assign(&state[n_cols + i]);
            state[i] = new_x;
        }
    }

    fn apply_mds_only(&self, state: &mut [F]) {
        let n_cols = self.params.n_cols;
        let x = self.mds_mul(&state[..n_cols]);
        let mut y = state[n_cols..].to_vec();
        y.rotate_left(1);
        let y = self.mds_mul(&y);

        state[..n_cols].clone_from_slice(&x);
        state[n_cols..].clone_from_slice(&y);
    }

    fn sbox_layer(&self, state: &mut [F]) {
        let n_cols = self.params.n_cols;
        for i in 0..n_cols {
            let mut x = state[i].clone();
            let mut y = state[n_cols + i].clone();

            let mut y_pow = y.clone();
            y_pow.square();
            let mut beta_y_pow = self.params.beta.clone();
            beta_y_pow.mul_assign(&y_pow);
            x.sub_assign(&beta_y_pow);

            let x_alpha_inv = x.pow_words_le(&self.params.alpha_inv);
            y.sub_assign(&x_alpha_inv);

            let mut y_pow_new = y.clone();
            y_pow_new.square();
            let mut beta_y_pow = self.params.beta.clone();
            beta_y_pow.mul_assign(&y_pow_new);
            x.add_assign(&beta_y_pow);
            x.add_assign(&self.params.delta);

            state[i] = x;
            state[n_cols + i] = y;
        }
    }
}

#[cfg(test)]
mod paired_profile_tests {
    use super::*;
    use crate::anemoi::instances::*;

    fn assert_same_non_round_profile<F: FieldElement>(
        current: &AnemoiParams<F>,
        original: &AnemoiParams<F>,
    ) {
        assert_eq!(current.n_cols, original.n_cols);
        assert_eq!(current.width, original.width);
        assert_eq!(current.alpha, original.alpha);
        assert_eq!(current.alpha_inv, original.alpha_inv);
        assert_eq!(current.beta, original.beta);
        assert_eq!(current.delta, original.delta);
        assert_eq!(current.mds, original.mds);
        assert_eq!(
            &current.round_constants_c[..original.rounds],
            original.round_constants_c
        );
        assert_eq!(
            &current.round_constants_d[..original.rounds],
            original.round_constants_d
        );
    }

    #[test]
    fn paired_profiles_share_all_non_round_choices() {
        assert_same_non_round_profile(
            &ANEMOI_BN254_4_CURRENT_PARAMS,
            &ANEMOI_BN254_4_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &ANEMOI_BLS12_381_4_CURRENT_PARAMS,
            &ANEMOI_BLS12_381_4_ORIGINAL_PARAMS,
        );
    }
}
