use crate::fields::FieldElement;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct ArionParams<F: FieldElement> {
    pub(crate) t: usize,
    /// High-degree exponent for the inverse S-box on the last state element.
    pub(crate) d: u64,
    /// d^{-1} mod (p-1).
    pub(crate) d_inv: [u64; 4],
    pub(crate) rounds: usize,
    /// GTDS g-coefficients: [round][non_last_idx][0] = alpha1, [1] = alpha2.
    /// g(s) = s^2 + alpha1 * s + alpha2 must be irreducible over F_p.
    pub(crate) g_coeffs: Vec<Vec<[F; 2]>>,
    /// GTDS h-coefficients: [round][non_last_idx] = beta1.
    /// h(s) = s^2 + beta1 * s.
    pub(crate) h_coeffs: Vec<Vec<F>>,
    /// Affine round constants applied after circulant matrix: [round][elem].
    pub(crate) round_constants: Vec<Vec<F>>,
}

impl<F: FieldElement> ArionParams<F> {
    pub fn new(
        t: usize,
        d: u64,
        d_inv: [u64; 4],
        rounds: usize,
        g_coeffs: &[Vec<[F; 2]>],
        h_coeffs: &[Vec<F>],
        round_constants: &[Vec<F>],
    ) -> Self {
        assert!(t >= 3);
        assert_eq!(g_coeffs.len(), rounds);
        assert_eq!(h_coeffs.len(), rounds);
        assert_eq!(round_constants.len(), rounds);
        for r in 0..rounds {
            assert_eq!(g_coeffs[r].len(), t - 1);
            assert_eq!(h_coeffs[r].len(), t - 1);
            assert_eq!(round_constants[r].len(), t);
        }
        assert_eq!(d, 257, "the canonical Arion high-degree exponent is 257");
        ArionParams {
            t,
            d,
            d_inv,
            rounds,
            g_coeffs: g_coeffs.to_owned(),
            h_coeffs: h_coeffs.to_owned(),
            round_constants: round_constants.to_owned(),
        }
    }

    pub fn exponent(&self) -> u64 {
        self.d
    }
}

#[cfg(test)]
mod tests {
    use super::{Arion, ArionParams};
    use crate::arion::instances::{ARION_BLS12_381_3_PARAMS, ARION_BN254_3_PARAMS};
    use crate::fields::bls12_381::Bls12_381;
    use crate::fields::bn254::Bn254;
    use crate::fields::FieldElement;
    use std::sync::Arc;

    fn fingerprint<F: FieldElement>(values: &[F]) -> F {
        let multiplier = F::from_u64(257);
        values.iter().fold(F::zero(), |mut acc, value| {
            acc.mul_assign(&multiplier);
            acc.add_assign(value);
            acc
        })
    }

    fn check_profile<F: FieldElement>(params: &Arc<ArionParams<F>>, expected_fingerprint: F) {
        assert_eq!(params.t, 3);
        assert_eq!(params.exponent(), 257);
        assert_eq!(params.rounds, 11);
        assert_eq!(params.g_coeffs.len(), 11);
        assert_eq!(params.h_coeffs.len(), 11);
        assert_eq!(params.round_constants.len(), 11);
        let input = (1..=3).map(F::from_u64).collect::<Vec<_>>();
        let output = Arion::new(params).permutation(&input);
        assert_eq!(fingerprint(&output), expected_fingerprint);
    }

    #[test]
    fn source_selected_t3_profiles_and_kats_are_exact() {
        check_profile(
            &ARION_BLS12_381_3_PARAMS,
            Bls12_381::from_hex("c7cfc02963665d0b576669d173b4348936ced9bc07d9ba12d13b04397668d1a")
                .unwrap(),
        );
        check_profile(
            &ARION_BN254_3_PARAMS,
            Bn254::from_hex("a076d031571bf3e56245be368a266dbf38e99c7e7fc7f2e78cb24b9790437b9")
                .unwrap(),
        );
    }
}

#[derive(Clone, Debug)]
pub struct Arion<F: FieldElement> {
    pub(crate) params: Arc<ArionParams<F>>,
}

impl<F: FieldElement> Arion<F> {
    pub fn new(params: &Arc<ArionParams<F>>) -> Self {
        Arion {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    // Arion permutation: (L_circ + c_R) ∘ GTDS ∘ ... ∘ (L_circ + c_0) ∘ GTDS.
    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut state = input.to_vec();

        for r in 0..self.params.rounds {
            self.gtds(&mut state, r);
            state = self.affine_layer(&state, r);
        }

        state
    }
    fn gtds(&self, state: &mut [F], round: usize) {
        let t = self.params.t;
        let mut f = vec![F::zero(); t];

        f[t - 1] = state[t - 1].pow_words_le(&self.params.d_inv);

        let mut sigma = state[t - 1].clone();
        sigma.add_assign(&f[t - 1]);

        for k in (0..t - 1).rev() {
            let x5 = state[k].pow_u64(5);

            // sigma^2
            let mut sigma_sq = sigma.clone();
            sigma_sq.square();

            // g(sigma) = sigma^2 + alpha1 * sigma + alpha2
            let g_coeff = &self.params.g_coeffs[round][k];
            let mut g = sigma_sq.clone();
            let mut alpha1_s = sigma.clone();
            alpha1_s.mul_assign(&g_coeff[0]);
            g.add_assign(&alpha1_s);
            g.add_assign(&g_coeff[1]);

            // h(sigma) = sigma^2 + beta1 * sigma
            let h_coeff = &self.params.h_coeffs[round][k];
            let mut h = sigma_sq;
            let mut beta1_s = sigma.clone();
            beta1_s.mul_assign(h_coeff);
            h.add_assign(&beta1_s);

            // f[k] = x^5 * g(sigma) + h(sigma)
            f[k] = x5;
            f[k].mul_assign(&g);
            f[k].add_assign(&h);

            // sigma += f[k] + state[k]
            sigma.add_assign(&f[k]);
            sigma.add_assign(&state[k]);
        }

        for (s, val) in state.iter_mut().zip(f.iter()) {
            *s = val.clone();
        }
    }

    fn apply_circulant(&self, x: &[F]) -> Vec<F> {
        let t = x.len();
        let sum = x.iter().fold(F::zero(), |mut acc, xi| {
            acc.add_assign(xi);
            acc
        });
        let mut out = vec![F::zero(); t];
        for i in 0..t {
            let mut acc = sum.clone();
            for (j, xj) in x.iter().enumerate() {
                let coeff = F::from_u64(((j + t - i) % t) as u64);
                let mut term = xj.clone();
                term.mul_assign(&coeff);
                acc.add_assign(&term);
            }
            out[i] = acc;
        }
        out
    }

    fn affine_layer(&self, x: &[F], round: usize) -> Vec<F> {
        let mut out = self.apply_circulant(x);
        let rc = &self.params.round_constants[round];
        for (o, c) in out.iter_mut().zip(rc.iter()) {
            o.add_assign(c);
        }
        out
    }
}
