use crate::fields::{biguint_to_limbs_le_4, PrimeFieldWords};
use crate::utils::read_field_from_shake;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake128;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct ReinforcedConcreteParams<F: PrimeFieldWords> {
    pub(crate) round_constants: Vec<Vec<F>>, // [round_idx][state_idx]
    pub(crate) alphas: Vec<u16>,
    pub(crate) betas: Vec<F>,
    pub(crate) si: Vec<u16>,
    pub(crate) sbox: Vec<u16>,
}

impl<F: PrimeFieldWords> ReinforcedConcreteParams<F> {
    pub const T: usize = 3;
    pub const PRE_ROUNDS: usize = 3;
    pub const POST_ROUNDS: usize = 3;
    pub const TOTAL_ROUNDS: usize = Self::PRE_ROUNDS + Self::POST_ROUNDS + 1;
    pub const INIT_SHAKE: &'static str = "ReinforcedConcrete";

    pub fn new(si: &[u16], sbox: &[u16], alphas: &[u16], betas: &[u16]) -> Self {
        assert_eq!(alphas.len(), Self::T - 1);
        assert_eq!(betas.len(), Self::T - 1);

        let mut shake = Self::init_shake();
        let round_constants = Self::instantiate_rc(&mut shake);

        let betas_field = betas.iter().map(|b| F::from_u64(*b as u64)).collect();

        ReinforcedConcreteParams {
            round_constants,
            alphas: alphas.to_owned(),
            betas: betas_field,
            si: si.to_owned(),
            sbox: Self::pad_sbox(sbox, si),
        }
    }

    // Domain-separate RC generation by algorithm tag and field modulus.
    fn init_shake() -> impl XofReader {
        let mut shake = Shake128::default();
        shake.update(Self::INIT_SHAKE.as_bytes());

        let limbs = biguint_to_limbs_le_4(&F::modulus());
        for limb in limbs {
            shake.update(&u64::to_le_bytes(limb));
        }

        shake.finalize_xof()
    }

    // Build the full round-constant matrix with shape [TOTAL_ROUNDS + 1][T].
    fn instantiate_rc(shake: &mut dyn XofReader) -> Vec<Vec<F>> {
        (0..=Self::TOTAL_ROUNDS)
            .map(|_| {
                (0..Self::T)
                    .map(|_| read_field_from_shake::<F>(shake))
                    .collect()
            })
            .collect()
    }

    // Pad the lookup table with identity values so Bars indexing is always defined.
    fn pad_sbox(sbox: &[u16], si: &[u16]) -> Vec<u16> {
        let max = *si.iter().max().expect("si must not be empty");
        assert!(sbox.len() <= max as usize);

        let mut out = sbox.to_owned();
        out.reserve((max as usize) - sbox.len());
        for i in (sbox.len() as u16)..max {
            out.push(i);
        }
        out
    }
}

#[derive(Clone, Debug)]
pub struct ReinforcedConcrete<F: PrimeFieldWords> {
    pub(crate) params: Arc<ReinforcedConcreteParams<F>>,
}

impl<F: PrimeFieldWords> ReinforcedConcrete<F> {
    pub fn new(params: &Arc<ReinforcedConcreteParams<F>>) -> Self {
        ReinforcedConcrete {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        ReinforcedConcreteParams::<F>::T
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        assert_eq!(input.len(), ReinforcedConcreteParams::<F>::T);

        let mut state = input.to_vec();
        self.concrete_in_place(&mut state, 0);

        for round in 1..=ReinforcedConcreteParams::<F>::PRE_ROUNDS {
            state = self.bricks(&state);
            self.concrete_in_place(&mut state, round);
        }

        state = self.bars(&state);
        self.concrete_in_place(&mut state, ReinforcedConcreteParams::<F>::PRE_ROUNDS + 1);

        for round in (ReinforcedConcreteParams::<F>::PRE_ROUNDS + 2)
            ..=ReinforcedConcreteParams::<F>::TOTAL_ROUNDS
        {
            state = self.bricks(&state);
            self.concrete_in_place(&mut state, round);
        }

        state
    }

    fn concrete_in_place(&self, state: &mut [F], round: usize) {
        // Multiplication by circulant(2,1,1,...) is state + sum(state).
        let mut sum = F::zero();
        for x in state.iter() {
            sum.add_assign(x);
        }

        for (x, rc) in state
            .iter_mut()
            .zip(self.params.round_constants[round].iter())
        {
            x.add_assign(&sum);
            x.add_assign(rc);
        }
    }

    fn bricks(&self, state: &[F]) -> Vec<F> {
        let mut out = state.to_vec();

        let mut x0_2 = state[0].clone();
        x0_2.square();
        x0_2.square();
        x0_2.mul_assign(&state[0]);
        out[0] = x0_2;

        for i in 1..state.len() {
            let prev = &state[i - 1];
            let mut prev_sq = prev.clone();
            prev_sq.square();

            for _ in 0..self.params.alphas[i - 1] {
                prev_sq.add_assign(prev);
            }
            prev_sq.add_assign(&self.params.betas[i - 1]);
            prev_sq.mul_assign(&state[i]);

            out[i] = prev_sq;
        }

        out
    }

    fn bars(&self, state: &[F]) -> Vec<F> {
        let mut out = state.to_vec();
        for x in out.iter_mut() {
            let mut digits = self.decompose(x);
            for d in digits.iter_mut() {
                *d = self.params.sbox[*d as usize];
            }
            *x = self.compose(&digits);
        }
        out
    }

    fn decompose(&self, val: &F) -> Vec<u16> {
        let len = self.params.si.len();
        let mut out = vec![0u16; len];
        let mut n = val.to_words_le();

        for i in (1..len).rev() {
            let (q, r) = div_rem_words_u16(n, self.params.si[i]);
            out[i] = r;
            n = q;
        }

        out[0] = words_to_u16_strict(n);

        out
    }

    fn compose(&self, vals: &[u16]) -> F {
        assert_eq!(vals.len(), self.params.si.len());

        let mut n = [0u64; 4];
        n[0] = vals[0] as u64;
        for (val, base) in vals.iter().zip(self.params.si.iter()).skip(1) {
            n = mul_add_words_u16(n, *base, *val);
        }

        F::from_words_le(n)
    }
}

fn div_rem_words_u16(mut n: [u64; 4], base: u16) -> ([u64; 4], u16) {
    assert!(base > 0, "division base must be non-zero");
    let base_u128 = base as u128;
    let mut rem = 0u128;

    for i in (0..4).rev() {
        let cur = (rem << 64) | n[i] as u128;
        n[i] = (cur / base_u128) as u64;
        rem = cur % base_u128;
    }

    (n, rem as u16)
}

fn mul_add_words_u16(n: [u64; 4], base: u16, add: u16) -> [u64; 4] {
    let base_u128 = base as u128;
    let mut out = [0u64; 4];
    let mut carry = add as u128;

    for i in 0..4 {
        let acc = (n[i] as u128) * base_u128 + carry;
        out[i] = acc as u64;
        carry = acc >> 64;
    }

    assert_eq!(carry, 0, "overflow while composing RC digits");
    out
}

fn words_to_u16_strict(n: [u64; 4]) -> u16 {
    assert_eq!(n[1], 0, "high limb must be zero in RC decomposition");
    assert_eq!(n[2], 0, "high limb must be zero in RC decomposition");
    assert_eq!(n[3], 0, "high limb must be zero in RC decomposition");
    u16::try_from(n[0]).expect("most significant RC digit must fit into u16")
}
