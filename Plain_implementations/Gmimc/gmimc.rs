use crate::fields::FieldElement;
use std::sync::Arc;

/// Parameters for one GMiMC expanding-round-function permutation.
///
/// Each round computes `y = (state[0] + c_r)^alpha`, adds `y` to every
/// other branch, and rotates the state left by one branch.
#[derive(Clone, Debug)]
pub struct GmimcParams<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) alpha: u64,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<F>,
}

impl<F: FieldElement> GmimcParams<F> {
    pub fn new(t: usize, alpha: u64, rounds: usize, round_constants: &[F]) -> Self {
        assert!(t > 1, "GMiMC needs at least two branches");
        assert!(alpha >= 2, "GMiMC alpha must be at least two");
        assert!(rounds > 0, "GMiMC needs at least one round");
        assert_eq!(
            round_constants.len(),
            rounds,
            "GMiMC needs exactly one constant per round"
        );
        Self {
            t,
            alpha,
            rounds,
            round_constants: round_constants.to_vec(),
        }
    }

    pub fn get_t(&self) -> usize {
        self.t
    }

    pub fn get_alpha(&self) -> u64 {
        self.alpha
    }

    pub fn get_rounds(&self) -> usize {
        self.rounds
    }
}

/// Plain GMiMC-erf permutation.
#[derive(Clone, Debug)]
pub struct Gmimc<F: FieldElement> {
    params: Arc<GmimcParams<F>>,
}

impl<F: FieldElement> Gmimc<F> {
    pub fn new(params: &Arc<GmimcParams<F>>) -> Self {
        Self {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn get_alpha(&self) -> u64 {
        self.params.alpha
    }

    pub fn get_rounds(&self) -> usize {
        self.params.rounds
    }

    #[inline(always)]
    fn sbox(&self, value: &F, round: usize) -> F {
        let mut x = value.clone();
        x.add_assign(&self.params.round_constants[round]);

        match self.params.alpha {
            3 => {
                let mut x2 = x.clone();
                x2.square();
                x2.mul_assign(&x);
                x2
            }
            5 => {
                let mut x4 = x.clone();
                x4.square();
                x4.square();
                x4.mul_assign(&x);
                x4
            }
            7 => {
                let mut x3 = x.clone();
                x3.square();
                x3.mul_assign(&x);
                let mut x7 = x3.clone();
                x7.square();
                x7.mul_assign(&x);
                x7
            }
            alpha => x.pow_u64(alpha),
        }
    }

    /// Apply the permutation and return a newly allocated state.
    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let mut state = input.to_vec();
        self.permutation_in_place(&mut state);
        state
    }

    /// Apply the accumulator implementation used for the wider small-prime states.
    pub fn permutation_in_place(&self, state: &mut [F]) {
        self.validate_state(state);
        if self.params.t < 8 {
            self.permutation_reference_in_place(state);
            return;
        }

        let t = self.params.t;
        let mut accumulator = F::zero();
        let mut window = vec![F::zero(); t - 1];

        for round in 0..self.params.rounds {
            let power = self.sbox(&state[0], round);

            window.rotate_left(1);
            accumulator.sub_assign(&window[0]);
            window[0] = power;
            accumulator.add_assign(&window[0]);

            state.rotate_left(1);
            state[0].add_assign(&accumulator);
        }

        // Settle the branches that are still in flight after the final round.
        for word in state.iter_mut().take(t - 1).skip(1) {
            window.rotate_left(1);
            accumulator.sub_assign(&window[0]);
            word.add_assign(&accumulator);
        }
    }

    /// Literal round-by-round oracle matching the sister repositories.
    pub fn permutation_reference(&self, input: &[F]) -> Vec<F> {
        let mut state = input.to_vec();
        self.validate_state(&state);
        self.permutation_reference_in_place(&mut state);
        state
    }

    fn permutation_reference_in_place(&self, state: &mut [F]) {
        for round in 0..self.params.rounds {
            let power = self.sbox(&state[0], round);
            for word in state.iter_mut().skip(1) {
                word.add_assign(&power);
            }
            state.rotate_left(1);
        }
    }

    fn validate_state(&self, state: &[F]) {
        assert_eq!(
            state.len(),
            self.params.t,
            "GMiMC input width does not match its parameters"
        );
    }
}
