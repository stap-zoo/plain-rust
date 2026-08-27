use super::instances::{LOOKUP_TABLE, NUM_SPLIT_AND_LOOKUP};
use crate::fields::goldilocks::Goldilocks;
use crate::fields::FieldElement;
use std::sync::Arc;

pub trait Tip5Field: FieldElement {
    fn to_u64(&self) -> u64;
}

impl Tip5Field for Goldilocks {
    fn to_u64(&self) -> u64 {
        Goldilocks::to_u64(self)
    }
}

#[derive(Clone, Debug)]
pub struct Tip5Params<F: Tip5Field> {
    pub(crate) t: usize,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<Vec<F>>,
    pub(crate) mds_first_column: Vec<F>,
    pub(crate) r: F,
    pub(crate) r_inv: F,
}

impl<F: Tip5Field> Tip5Params<F> {
    pub fn new(
        t: usize,
        rounds: usize,
        round_constants: Vec<Vec<F>>,
        mds_first_column: Vec<F>,
        r: F,
        r_inv: F,
    ) -> Self {
        Tip5Params {
            t,
            rounds,
            round_constants,
            mds_first_column,
            r,
            r_inv,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Tip5<F: Tip5Field> {
    pub(crate) params: Arc<Tip5Params<F>>,
}

impl<F: Tip5Field> Tip5<F> {
    pub fn new(params: &Arc<Tip5Params<F>>) -> Self {
        Tip5 {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut state = input.to_vec();
        for round in 0..self.params.rounds {
            self.sbox_layer(&mut state);
            state = self.mds_matmul(&state);
            self.add_round_constants(&mut state, round);
        }
        state
    }

    fn add_round_constants(&self, state: &mut [F], round: usize) {
        for (el, rc) in state
            .iter_mut()
            .zip(self.params.round_constants[round].iter())
        {
            el.add_assign(rc);
        }
    }

    fn sbox_layer(&self, state: &mut [F]) {
        for i in 0..NUM_SPLIT_AND_LOOKUP {
            state[i] = self.split_and_lookup(&state[i]);
        }
        for i in NUM_SPLIT_AND_LOOKUP..self.params.t {
            state[i] = state[i].pow_u64(7);
        }
    }

    fn split_and_lookup(&self, element: &F) -> F {
        let mut monty = element.clone();
        monty.mul_assign(&self.params.r);
        let mut bytes = monty.to_u64().to_le_bytes();
        for b in bytes.iter_mut() {
            *b = LOOKUP_TABLE[*b as usize];
        }
        let mut out = F::from_u64(u64::from_le_bytes(bytes));
        out.mul_assign(&self.params.r_inv);
        out
    }

    /// Circulant MDS multiplication using first-column indexing.
    /// M[row][col] = mds_first_column[(row - col + t) % t].
    /// Matches the original twenty-first circulant MDS implementation.
    fn mds_matmul(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        debug_assert_eq!(t, input.len());
        let col = &self.params.mds_first_column;
        let mut out = vec![F::zero(); t];
        for row in 0..t {
            for (col_idx, inp) in input.iter().enumerate() {
                let idx = (t + row - col_idx) % t;
                let mut tmp = col[idx].clone();
                tmp.mul_assign(inp);
                out[row].add_assign(&tmp);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::super::instances::TIP5_GOLDILOCKS_PARAMS;
    use super::Tip5;
    use crate::fields::goldilocks::Goldilocks;
    use crate::fields::FieldElement;

    #[test]
    fn permutation_matches_tip5_hash10_reference_vector_0() {
        let perm = Tip5::new(&TIP5_GOLDILOCKS_PARAMS);
        let mut input = vec![Goldilocks::from_u64(0); perm.get_t()];
        for x in input.iter_mut().skip(10) {
            *x = Goldilocks::from_u64(1);
        }

        let output = perm.permutation(&input);
        let expected = [
            941080798860502477u64,
            5295886365985465639u64,
            14728839126885177993u64,
            10358449902914633406u64,
            14220746792122877272u64,
        ];

        for (got, want) in output.iter().take(expected.len()).zip(expected.iter()) {
            assert_eq!(got.to_u64(), *want);
        }
    }
}
