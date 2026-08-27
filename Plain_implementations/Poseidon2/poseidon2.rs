use crate::fields::FieldElement;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct Poseidon2Params<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) d: u64,
    pub(crate) rounds_f_beginning: usize,
    pub(crate) rounds_p: usize,
    pub(crate) rounds: usize,
    pub(crate) mat_internal_diag_m_1: Vec<F>,
    pub(crate) mat_external: Vec<Vec<F>>,
    pub(crate) mat_internal: Vec<Vec<F>>,
    pub(crate) round_constants: Vec<Vec<F>>, // [round_idx][state_idx]
    // The IAIK BLS12-381 t=4 profile uses 2*M4, while the published
    // gnark-crypto/HorizenLabs BN254 t=4 profile uses M4.
    pub(crate) t4_external_is_m4: bool,
}

impl<F: FieldElement> Poseidon2Params<F> {
    pub fn new(
        t: usize,
        d: u64,
        rounds_f: usize,
        rounds_p: usize,
        mat_external: &[Vec<F>],
        mat_internal: &[Vec<F>],
        round_constants: &[Vec<F>],
    ) -> Self {
        Self::new_with_t4_external_convention(
            t,
            d,
            rounds_f,
            rounds_p,
            mat_external,
            mat_internal,
            round_constants,
            false,
        )
    }

    pub(crate) fn new_with_t4_m4_external(
        t: usize,
        d: u64,
        rounds_f: usize,
        rounds_p: usize,
        mat_external: &[Vec<F>],
        mat_internal: &[Vec<F>],
        round_constants: &[Vec<F>],
    ) -> Self {
        assert_eq!(t, 4, "the M4-only external convention is specific to t=4");
        Self::new_with_t4_external_convention(
            t,
            d,
            rounds_f,
            rounds_p,
            mat_external,
            mat_internal,
            round_constants,
            true,
        )
    }

    fn new_with_t4_external_convention(
        t: usize,
        d: u64,
        rounds_f: usize,
        rounds_p: usize,
        mat_external: &[Vec<F>],
        mat_internal: &[Vec<F>],
        round_constants: &[Vec<F>],
        t4_external_is_m4: bool,
    ) -> Self {
        let mut mat_internal_diag_m_1 = Vec::with_capacity(t);
        let one = F::one();
        for (i, row) in mat_internal.iter().enumerate() {
            let mut diag_m_1 = row[i].clone();
            diag_m_1.sub_assign(&one);
            mat_internal_diag_m_1.push(diag_m_1);
        }

        Poseidon2Params {
            t,
            d,
            rounds_f_beginning: rounds_f / 2,
            rounds_p,
            rounds: rounds_f + rounds_p,
            mat_internal_diag_m_1,
            mat_external: mat_external.to_owned(),
            mat_internal: mat_internal.to_owned(),
            round_constants: round_constants.to_owned(),
            t4_external_is_m4,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Poseidon2<F: FieldElement> {
    pub(crate) params: Arc<Poseidon2Params<F>>,
}

impl<F: FieldElement> Poseidon2<F> {
    pub fn new(params: &Arc<Poseidon2Params<F>>) -> Self {
        Poseidon2 {
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

        self.matmul_external(&mut state);

        for r in 0..self.params.rounds_f_beginning {
            self.add_rc_in_place(&mut state, r);
            self.sbox_in_place(&mut state);
            self.matmul_external(&mut state);
        }

        let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
        for r in self.params.rounds_f_beginning..p_end {
            state[0].add_assign(&self.params.round_constants[r][0]);
            state[0] = self.sbox_p(&state[0]);
            self.matmul_internal(&mut state);
        }

        for r in p_end..self.params.rounds {
            self.add_rc_in_place(&mut state, r);
            self.sbox_in_place(&mut state);
            self.matmul_external(&mut state);
        }

        state
    }

    fn sbox_in_place(&self, state: &mut [F]) {
        for el in state.iter_mut() {
            *el = self.sbox_p(el);
        }
    }

    fn sbox_p(&self, input: &F) -> F {
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

    fn matmul_external(&self, input: &mut [F]) {
        let t = self.params.t;
        match t {
            2 => {
                let mut sum = input[0].clone();
                sum.add_assign(&input[1]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
            }
            3 => {
                let mut sum = input[0].clone();
                sum.add_assign(&input[1]);
                sum.add_assign(&input[2]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
                input[2].add_assign(&sum);
            }
            4 | 8 | 12 | 16 | 24 => {
                let t4 = t / 4;
                for i in 0..t4 {
                    let start_index = i * 4;
                    let mut t_0 = input[start_index].clone();
                    t_0.add_assign(&input[start_index + 1]);
                    let mut t_1 = input[start_index + 2].clone();
                    t_1.add_assign(&input[start_index + 3]);
                    let mut t_2 = input[start_index + 1].clone();
                    t_2.double();
                    t_2.add_assign(&t_1);
                    let mut t_3 = input[start_index + 3].clone();
                    t_3.double();
                    t_3.add_assign(&t_0);
                    let mut t_4 = t_1.clone();
                    t_4.double();
                    t_4.double();
                    t_4.add_assign(&t_3);
                    let mut t_5 = t_0.clone();
                    t_5.double();
                    t_5.double();
                    t_5.add_assign(&t_2);
                    let mut t_6 = t_3.clone();
                    t_6.add_assign(&t_5);
                    let mut t_7 = t_2.clone();
                    t_7.add_assign(&t_4);
                    input[start_index] = t_6;
                    input[start_index + 1] = t_5;
                    input[start_index + 2] = t_7;
                    input[start_index + 3] = t_4;
                }

                if t == 4 && self.params.t4_external_is_m4 {
                    return;
                }

                let mut stored = vec![F::zero(); 4];
                for l in 0..4 {
                    stored[l] = input[l].clone();
                    for j in 1..t4 {
                        stored[l].add_assign(&input[4 * j + l]);
                    }
                }
                for i in 0..input.len() {
                    input[i].add_assign(&stored[i % 4]);
                }
            }
            _ => panic!("unsupported width"),
        }
    }

    fn matmul_internal(&self, input: &mut [F]) {
        let t = self.params.t;

        match t {
            2 => {
                let mut sum = input[0].clone();
                sum.add_assign(&input[1]);
                input[0].add_assign(&sum);
                input[1].double();
                input[1].add_assign(&sum);
            }
            3 => {
                let mut sum = input[0].clone();
                sum.add_assign(&input[1]);
                sum.add_assign(&input[2]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
                input[2].double();
                input[2].add_assign(&sum);
            }
            4 | 8 | 12 | 16 | 24 => {
                let mut sum = input[0].clone();
                for el in input.iter().skip(1) {
                    sum.add_assign(el);
                }
                for i in 0..input.len() {
                    input[i].mul_assign(&self.params.mat_internal_diag_m_1[i]);
                    input[i].add_assign(&sum);
                }
            }
            _ => panic!("unsupported width"),
        }
    }

    fn add_rc_in_place(&self, state: &mut [F], round: usize) {
        let rc = &self.params.round_constants[round];
        for (x, c) in state.iter_mut().zip(rc.iter()) {
            x.add_assign(c);
        }
    }
}
