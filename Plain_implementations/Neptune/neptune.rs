use crate::fields::FieldElement;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct NeptuneParams<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) d: u64,
    pub(crate) rounds_f_beginning: usize,
    pub(crate) rounds_p: usize,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<Vec<F>>, // [round_idx][state_idx]
    pub(crate) m_e: Vec<Vec<F>>,             // external matrix
    pub(crate) mu: Vec<F>,                   // diagonal minus one of internal matrix
    pub(crate) gamma: F,
}

impl<F: FieldElement> NeptuneParams<F> {
    pub fn new(
        t: usize,
        d: u64,
        rounds_f: usize,
        rounds_p: usize,
        m_e: &[Vec<F>],
        mu: &[F],
        round_constants: &[Vec<F>],
        gamma: &F,
    ) -> Self {
        assert!(matches!(d, 3 | 5 | 7));
        assert_eq!(rounds_f % 2, 0);
        assert_eq!(t % 2, 0);
        assert_eq!(m_e.len(), t);
        assert_eq!(mu.len(), t);
        assert_eq!(round_constants.len(), rounds_f + rounds_p);
        for row in m_e {
            assert_eq!(row.len(), t);
        }
        for rc in round_constants {
            assert_eq!(rc.len(), t);
        }

        NeptuneParams {
            t,
            d,
            rounds_f_beginning: rounds_f / 2,
            rounds_p,
            rounds: rounds_f + rounds_p,
            round_constants: round_constants.to_owned(),
            m_e: m_e.to_owned(),
            mu: mu.to_owned(),
            gamma: gamma.clone(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Neptune<F: FieldElement> {
    pub(crate) params: Arc<NeptuneParams<F>>,
}

impl<F: FieldElement> Neptune<F> {
    pub fn new(params: &Arc<NeptuneParams<F>>) -> Self {
        Neptune {
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
        self.external_matmul(&mut state);

        let half_f = self.params.rounds_f_beginning;
        let mut round = 0usize;

        for _ in 0..half_f {
            self.external_round(&mut state, round);
            round += 1;
        }

        for _ in 0..self.params.rounds_p {
            self.internal_round(&mut state, round);
            round += 1;
        }

        while round < self.params.rounds {
            self.external_round(&mut state, round);
            round += 1;
        }

        state
    }

    #[inline(always)]
    fn external_round(&self, state: &mut [F], round: usize) {
        self.external_sbox_in_place(state);
        self.external_matmul(state);
        self.add_rc_in_place(state, round);
    }

    #[inline(always)]
    fn internal_round(&self, state: &mut [F], round: usize) {
        state[0] = self.sbox_d(&state[0]);
        self.internal_matmul(state);
        self.add_rc_in_place(state, round);
    }

    fn sbox_d(&self, input: &F) -> F {
        let mut input2 = input.clone();
        input2.square();

        match self.params.d {
            3 => {
                let mut out = input2;
                out.mul_assign(input);
                out
            }
            5 => {
                let mut out = input2;
                out.square();
                out.mul_assign(input);
                out
            }
            7 => {
                let mut out = input2.clone();
                out.square();
                out.mul_assign(&input2);
                out.mul_assign(input);
                out
            }
            _ => panic!("unsupported s-box degree"),
        }
    }

    fn external_sbox_in_place(&self, state: &mut [F]) {
        for i in (0..state.len()).step_by(2) {
            let (y1, y2) = self.external_sbox_prime(&state[i], &state[i + 1]);
            state[i] = y1;
            state[i + 1] = y2;
        }
    }

    fn external_sbox_prime(&self, x1: &F, x2: &F) -> (F, F) {
        let mut zi = x1.clone();
        zi.sub_assign(x2);

        let mut zib = zi.clone();
        zib.square();

        let mut sum = x1.clone();
        sum.add_assign(x2);

        let mut y1 = sum.clone();
        y1.add_assign(x1);
        let mut y2 = sum.clone();
        y2.add_assign(x2);
        y2.add_assign(x2);

        let mut tmp1 = zib.clone();
        tmp1.double();
        let mut tmp2 = tmp1.clone();
        tmp1.add_assign(&zib);
        tmp2.double();
        y1.add_assign(&tmp1);
        y2.add_assign(&tmp2);

        let mut tmp = zi;
        tmp.sub_assign(x2);
        tmp.sub_assign(&zib);
        tmp.add_assign(&self.params.gamma);
        tmp.square();
        y1.add_assign(&tmp);
        y2.add_assign(&tmp);

        (y1, y2)
    }

    fn external_matmul(&self, state: &mut [F]) {
        match self.params.t {
            4 => {
                let out = Self::external_matmul_4(state);
                state.clone_from_slice(&out);
            }
            8 => {
                let out = Self::external_matmul_8(state);
                state.clone_from_slice(&out);
            }
            _ => {
                self.matmul_in_place(state, &self.params.m_e);
            }
        }
    }

    /// Hand-optimized circulant expansion for t=4: circ(2, 1, 1, 1).
    /// From zkfriendlyhashzoo: swap + sum + accumulate.
    fn external_matmul_4(input: &[F]) -> Vec<F> {
        let mut output = input.to_vec();
        output.swap(1, 3);

        let mut sum1 = input[0].clone();
        sum1.add_assign(&input[2]);
        let mut sum2 = input[1].clone();
        sum2.add_assign(&input[3]);

        output[0].add_assign(&sum1);
        output[1].add_assign(&sum2);
        output[2].add_assign(&sum1);
        output[3].add_assign(&sum2);

        output
    }

    /// Hand-optimized circulant expansion for t=8: circ(3, 2, 1, 1, ...).
    /// Equals state + state + rot(state) + sum(state).
    /// From zkfriendlyhashzoo: swap + sum + rotate + double.
    fn external_matmul_8(input: &[F]) -> Vec<F> {
        let mut output = input.to_vec();
        output.swap(1, 7);
        output.swap(3, 5);

        let mut sum1 = input[0].clone();
        let mut sum2 = input[1].clone();

        for el in input.iter().step_by(2).skip(1) {
            sum1.add_assign(el);
        }
        for el in input.iter().skip(1).step_by(2).skip(1) {
            sum2.add_assign(el);
        }

        let mut output_rot = output.clone();
        output_rot.rotate_left(2);

        for (i, el) in output.iter_mut().enumerate() {
            el.double();
            el.add_assign(&output_rot[i]);
            if i & 1 == 0 {
                el.add_assign(&sum1);
            } else {
                el.add_assign(&sum2);
            }
        }

        output.swap(3, 7);
        output
    }

    fn internal_matmul(&self, state: &mut [F]) {
        let mut sum = F::zero();
        for x in state.iter() {
            sum.add_assign(x);
        }

        for (x, mu) in state.iter_mut().zip(self.params.mu.iter()) {
            let mut out = x.clone();
            out.mul_assign(mu);
            out.add_assign(&sum);
            *x = out;
        }
    }

    fn matmul_in_place(&self, state: &mut [F], mat: &[Vec<F>]) {
        let t = state.len();
        let mut out = vec![F::zero(); t];
        for row in 0..t {
            let mut acc = F::zero();
            for (col, val) in state.iter().enumerate() {
                let mut tmp = mat[row][col].clone();
                tmp.mul_assign(val);
                acc.add_assign(&tmp);
            }
            out[row] = acc;
        }
        state.clone_from_slice(&out);
    }

    fn add_rc_in_place(&self, state: &mut [F], round: usize) {
        let rc = &self.params.round_constants[round];
        for (x, c) in state.iter_mut().zip(rc.iter()) {
            x.add_assign(c);
        }
    }
}
