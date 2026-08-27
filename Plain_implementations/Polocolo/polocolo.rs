use crate::fields::{FieldElement, PrimeFieldWords};
use std::collections::HashMap;
use std::sync::{Arc, OnceLock};

#[derive(Clone, Debug)]
pub struct PolocoloParams<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) m: usize,
    pub(crate) ann: [u64; 4],
    pub(crate) rounds: usize,
    pub(crate) mds: Vec<Vec<F>>,
    pub(crate) round_constants: Vec<Vec<F>>,
    pub(crate) lut: HashMap<[u64; 4], F>,
    pub(crate) modulus_minus_two: [u64; 4],
}

impl<F: FieldElement> PolocoloParams<F> {
    pub fn new(
        t: usize,
        m: usize,
        ann: [u64; 4],
        rounds: usize,
        mds: &[Vec<F>],
        round_constants: &[Vec<F>],
        lut: &HashMap<[u64; 4], F>,
        modulus_minus_two: [u64; 4],
    ) -> Self {
        assert_eq!(mds.len(), t);
        assert_eq!(round_constants.len(), rounds);
        assert_eq!(round_constants[0].len(), t);

        PolocoloParams {
            t,
            m,
            ann,
            rounds,
            mds: mds.to_owned(),
            round_constants: round_constants.to_owned(),
            lut: lut.to_owned(),
            modulus_minus_two,
        }
    }
}

#[derive(Debug)]
pub struct Polocolo<F: FieldElement> {
    pub(crate) params: Arc<PolocoloParams<F>>,
    fallback_lut: OnceLock<HashMap<[u64; 4], F>>,
    fallback_poly: OnceLock<Vec<F>>,
}

impl<F: FieldElement> Clone for Polocolo<F> {
    fn clone(&self) -> Self {
        Self {
            params: Arc::clone(&self.params),
            fallback_lut: OnceLock::new(),
            fallback_poly: OnceLock::new(),
        }
    }
}

impl<F: PrimeFieldWords> Polocolo<F> {
    pub fn new(params: &Arc<PolocoloParams<F>>) -> Self {
        Polocolo {
            params: Arc::clone(params),
            fallback_lut: OnceLock::new(),
            fallback_poly: OnceLock::new(),
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
            state = self.affine(&state, round);
            state = self.sbox(&state);
        }
        self.affine(&state, self.params.rounds)
    }

    fn sbox(&self, input: &[F]) -> Vec<F> {
        input.iter().map(|el| self.sbox_elem(el)).collect()
    }

    fn sbox_elem(&self, el: &F) -> F {
        if *el == F::zero() {
            return F::zero();
        }

        let lut_in = el.pow_words_le(&self.params.ann);
        let lut_key = lut_in.to_words_le();
        let lut_out = if let Some(out) = self.params.lut.get(&lut_key) {
            out
        } else if let Some(out) = self.recovered_lut().get(&lut_key) {
            out
        } else {
            let coeffs = self.fallback_poly();
            return self.sbox_with_poly(el, coeffs);
        };

        let mut out = el.pow_words_le(&self.params.modulus_minus_two);
        out.mul_assign(lut_out);
        out
    }

    fn recovered_lut(&self) -> &HashMap<[u64; 4], F> {
        self.fallback_lut
            .get_or_init(|| recover_missing_lut_entries(self.params.as_ref()))
    }

    fn fallback_poly(&self) -> &Vec<F> {
        self.fallback_poly
            .get_or_init(|| build_lut_interpolant(self.params.as_ref()))
    }

    fn sbox_with_poly(&self, el: &F, coeffs: &[F]) -> F {
        let lut_in = el.pow_words_le(&self.params.ann);
        let lut_out = eval_poly(coeffs, &lut_in);
        let mut out = el.pow_words_le(&self.params.modulus_minus_two);
        out.mul_assign(&lut_out);
        out
    }

    fn affine(&self, input: &[F], round: usize) -> Vec<F> {
        let mat_result = self.matmul(input, &self.params.mds);
        if round < self.params.rounds {
            Self::add_rc(&mat_result, &self.params.round_constants[round])
        } else {
            mat_result
        }
    }

    fn matmul(&self, input: &[F], mat: &[Vec<F>]) -> Vec<F> {
        let t = mat.len();
        debug_assert_eq!(t, input.len());

        let mut out = vec![F::zero(); t];
        for row in 0..t {
            for (col, inp) in input.iter().enumerate().take(t) {
                let mut tmp = mat[row][col].clone();
                tmp.mul_assign(inp);
                out[row].add_assign(&tmp);
            }
        }
        out
    }

    fn add_rc(input: &[F], round_constants: &[F]) -> Vec<F> {
        debug_assert_eq!(input.len(), round_constants.len());
        input
            .iter()
            .zip(round_constants.iter())
            .map(|(a, b)| {
                let mut r = a.clone();
                r.add_assign(b);
                r
            })
            .collect()
    }
}

fn recover_missing_lut_entries<F: PrimeFieldWords>(
    params: &PolocoloParams<F>,
) -> HashMap<[u64; 4], F> {
    let mut recovered = HashMap::new();
    let m = params.m;
    if m == 0 {
        return recovered;
    }

    let g = F::from_biguint(&F::generator());
    let gk = g.pow_words_le(&params.ann);

    let mut subgroup_index = HashMap::with_capacity(m);
    let mut h = F::one();
    for r in 0..m {
        subgroup_index.insert(h.to_words_le(), r);
        h.mul_assign(&gk);
    }

    let mut g_pow_index = HashMap::with_capacity(m);
    let mut gs = F::one();
    for s in 0..m {
        g_pow_index.insert(gs.to_words_le(), s);
        gs.mul_assign(&g);
    }

    let mut sigma: Vec<Option<usize>> = vec![None; m];
    let zero_key = F::zero().to_words_le();

    for (key, out) in params.lut.iter() {
        if *key == zero_key {
            continue;
        }
        let Some(&r) = subgroup_index.get(key) else {
            continue;
        };

        let base_exp = (m as u64 + 1) * (r as u64);
        let base = g.pow_u64(base_exp);
        let inv_base = base.pow_words_le(&params.modulus_minus_two);
        let mut ratio = out.clone();
        ratio.mul_assign(&inv_base);

        if let Some(&s) = g_pow_index.get(&ratio.to_words_le()) {
            sigma[r] = Some(s);
        }
    }

    let mut used_sigma = vec![false; m];
    let mut missing_r = Vec::new();
    for (r, maybe_sigma) in sigma.iter().enumerate() {
        if let Some(s) = maybe_sigma {
            used_sigma[*s] = true;
        } else {
            missing_r.push(r);
        }
    }

    let mut missing_sigma = Vec::new();
    for (s, used) in used_sigma.iter().enumerate() {
        if !used {
            missing_sigma.push(s);
        }
    }

    if missing_r.len() == 1 && missing_sigma.len() == 1 {
        sigma[missing_r[0]] = Some(missing_sigma[0]);
    } else if missing_r.len() == missing_sigma.len() {
        for (r, s) in missing_r.iter().zip(missing_sigma.iter()) {
            sigma[*r] = Some(*s);
        }
    } else {
        return recovered;
    }

    for r in missing_r {
        let Some(s) = sigma[r] else {
            continue;
        };

        let key = gk.pow_u64(r as u64).to_words_le();
        if params.lut.contains_key(&key) {
            continue;
        }

        let out_exp = (m as u64 + 1) * (r as u64) + s as u64;
        recovered.insert(key, g.pow_u64(out_exp));
    }

    recovered
}

fn build_lut_interpolant<F: PrimeFieldWords>(params: &PolocoloParams<F>) -> Vec<F> {
    let mut xs = Vec::with_capacity(params.lut.len());
    let mut ys = Vec::with_capacity(params.lut.len());
    for (x_words, y) in params.lut.iter() {
        xs.push(F::from_words_le(*x_words));
        ys.push(y.clone());
    }
    interpolate_poly(&xs, &ys, &params.modulus_minus_two)
}

fn eval_poly<F: FieldElement>(coeffs: &[F], x: &F) -> F {
    if coeffs.is_empty() {
        return F::zero();
    }
    let mut acc = coeffs[coeffs.len() - 1].clone();
    for coeff in coeffs[..coeffs.len() - 1].iter().rev() {
        acc.mul_assign(x);
        acc.add_assign(coeff);
    }
    acc
}

fn interpolate_poly<F: FieldElement>(xs: &[F], ys: &[F], modulus_minus_two: &[u64; 4]) -> Vec<F> {
    let m = xs.len();
    assert_eq!(ys.len(), m);

    let mut poly = vec![F::one()];
    for x in xs.iter() {
        poly = poly_mul_monic_linear(&poly, x);
    }

    let mut coeffs = vec![F::zero(); m];
    for i in 0..m {
        let mut denom = F::one();
        for j in 0..m {
            if i == j {
                continue;
            }
            let mut diff = xs[i].clone();
            diff.sub_assign(&xs[j]);
            denom.mul_assign(&diff);
        }
        let denom_inv = denom.pow_words_le(modulus_minus_two);
        let q = poly_div_monic_linear(&poly, &xs[i]);

        let mut scale = ys[i].clone();
        scale.mul_assign(&denom_inv);
        for k in 0..m {
            let mut term = q[k].clone();
            term.mul_assign(&scale);
            coeffs[k].add_assign(&term);
        }
    }

    coeffs
}

fn poly_mul_monic_linear<F: FieldElement>(poly: &[F], root: &F) -> Vec<F> {
    let mut out = vec![F::zero(); poly.len() + 1];
    for (i, coeff) in poly.iter().enumerate() {
        let mut neg = root.clone();
        neg.mul_assign(coeff);
        out[i].sub_assign(&neg);
        out[i + 1].add_assign(coeff);
    }
    out
}

fn poly_div_monic_linear<F: FieldElement>(poly: &[F], root: &F) -> Vec<F> {
    let m = poly.len() - 1;
    let mut q = vec![F::zero(); m];
    q[m - 1] = poly[m].clone();
    for k in (1..m).rev() {
        let mut term = root.clone();
        term.mul_assign(&q[k]);
        let mut val = poly[k].clone();
        val.add_assign(&term);
        q[k - 1] = val;
    }
    q
}
