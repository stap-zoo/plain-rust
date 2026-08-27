use super::instances::{LOOKUP_TABLE, NUM_SPLIT_AND_LOOKUP};
use crate::fields::goldilocks::Goldilocks;
use crate::fields::FieldElement;
use std::sync::Arc;

pub trait Tip4Field: FieldElement {
    fn to_u64(&self) -> u64;
}

impl Tip4Field for Goldilocks {
    fn to_u64(&self) -> u64 {
        Goldilocks::to_u64(self)
    }
}

#[derive(Clone, Debug)]
pub struct Tip4Params<F: Tip4Field> {
    pub(crate) t: usize,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<Vec<F>>,
    pub(crate) mds_first_column: Vec<F>,
    pub(crate) r: F,
    pub(crate) r_inv: F,
}

impl<F: Tip4Field> Tip4Params<F> {
    pub fn new(
        t: usize,
        rounds: usize,
        round_constants: Vec<Vec<F>>,
        mds_first_column: Vec<F>,
        r: F,
        r_inv: F,
    ) -> Self {
        Tip4Params {
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
pub struct Tip4<F: Tip4Field> {
    pub(crate) params: Arc<Tip4Params<F>>,
}

impl<F: Tip4Field> Tip4<F> {
    pub fn new(params: &Arc<Tip4Params<F>>) -> Self {
        Tip4 {
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
    use super::super::instances::TIP4P_GOLDILOCKS_PARAMS;
    use super::Tip4;
    use crate::fields::goldilocks::Goldilocks;
    use crate::fields::FieldElement;

    #[test]
    fn permutation_matches_tip4p_winterfell_reference_vector() {
        let perm = Tip4::new(&TIP4P_GOLDILOCKS_PARAMS);
        let input: Vec<Goldilocks> = (0..perm.get_t())
            .map(|i| Goldilocks::from_u64(i as u64))
            .collect();
        let output = perm.permutation(&input);

        let expected = [
            8086224146445274039u64,
            12620228105612859910u64,
            4429745645163147655u64,
            12827206290147492018u64,
            7103575185686863209u64,
            5938996934280238338u64,
            7458235737397060283u64,
            127950926479970750u64,
            433935175963827303u64,
            11405496933068372192u64,
            4026696970861104429u64,
            6779880475047698803u64,
        ];

        assert_eq!(output.len(), expected.len());
        for (got, want) in output.iter().zip(expected.iter()) {
            assert_eq!(got.to_u64(), *want);
        }
    }
}
