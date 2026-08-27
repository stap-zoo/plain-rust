use crate::fields::FieldElement;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct GriffinParams<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) d: u64,
    pub(crate) d_inv: [u64; 4],
    pub(crate) rounds: usize,
    pub(crate) alpha_beta: Vec<[F; 2]>,
    pub(crate) round_constants: Vec<Vec<F>>, // [round_idx][state_idx], only for first rounds-1 rounds
}

impl<F: FieldElement> GriffinParams<F> {
    pub fn new(
        t: usize,
        d: u64,
        d_inv: [u64; 4],
        rounds: usize,
        alpha_beta: &[[F; 2]],
        round_constants: &[Vec<F>],
    ) -> Self {
        assert!(t == 3 || (t >= 4 && t % 4 == 0));
        assert!(rounds >= 1);
        assert_eq!(alpha_beta.len(), t.saturating_sub(2));
        assert_eq!(round_constants.len(), rounds.saturating_sub(1));
        for rc in round_constants {
            assert_eq!(rc.len(), t);
        }

        GriffinParams {
            t,
            d,
            d_inv,
            rounds,
            alpha_beta: alpha_beta.to_vec(),
            round_constants: round_constants.to_owned(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Griffin<F: FieldElement> {
    pub(crate) params: Arc<GriffinParams<F>>,
}

impl<F: FieldElement> Griffin<F> {
    pub fn new(params: &Arc<GriffinParams<F>>) -> Self {
        Griffin {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn get_alpha(&self) -> u64 {
        self.params.d
    }

    pub fn get_rounds(&self) -> usize {
        self.params.rounds
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut state = input.to_vec();
        // Griffin-π: Gπ(x) = F_{R-1} ∘ ... ∘ F0(M * x)
        self.linear_layer(&mut state);

        for round in 0..self.params.rounds {
            state = self.non_linear(&state);
            self.linear_layer(&mut state);
            if round + 1 < self.params.rounds {
                self.add_rc_in_place(&mut state, round);
            }
        }

        state
    }

    fn non_linear(&self, input: &[F]) -> Vec<F> {
        let mut output = input.to_vec();

        output[0] = output[0].pow_words_le(&self.params.d_inv);
        output[1] = self.sbox_d(&input[1]);

        let y0 = output[0].clone();
        let mut y01 = y0.clone();
        y01.add_assign(&output[1]);

        for (i, ((out, prev), ab)) in output
            .iter_mut()
            .skip(2)
            .zip(input.iter().skip(1))
            .zip(self.params.alpha_beta.iter())
            .enumerate()
        {
            let l = if i == 0 {
                y01.clone()
            } else {
                y01.add_assign(&y0);
                let mut tmp = y01.clone();
                tmp.add_assign(prev);
                tmp
            };

            let mut poly = l.clone();
            poly.square();

            let mut alpha_l = l;
            alpha_l.mul_assign(&ab[0]);
            poly.add_assign(&alpha_l);
            poly.add_assign(&ab[1]);

            out.mul_assign(&poly);
        }

        output
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
            _ => input.pow_u64(self.params.d),
        }
    }

    fn linear_layer(&self, state: &mut [F]) {
        match state.len() {
            3 => {
                // circulant(2,1,1): add global sum to each coordinate.
                let mut sum = state[0].clone();
                for x in state.iter().skip(1) {
                    sum.add_assign(x);
                }
                for x in state.iter_mut() {
                    x.add_assign(&sum);
                }
            }
            4 => Self::apply_m4(state),
            t if t >= 8 && t % 4 == 0 => {
                // Structured multiplication by M = M' * M'' for t = 4t' >= 8.
                let t4 = t / 4;

                for chunk in state.chunks_exact_mut(4) {
                    Self::apply_m4(chunk);
                }

                let mut lanes = [F::zero(), F::zero(), F::zero(), F::zero()];
                for lane in 0..4 {
                    lanes[lane] = state[lane].clone();
                    for block in 1..t4 {
                        lanes[lane].add_assign(&state[4 * block + lane]);
                    }
                }

                for i in 0..state.len() {
                    state[i].add_assign(&lanes[i % 4]);
                }
            }
            _ => panic!("unsupported width"),
        }
    }

    #[inline(always)]
    fn apply_m4(state: &mut [F]) {
        debug_assert_eq!(state.len(), 4);

        let mut t0 = state[0].clone();
        t0.add_assign(&state[1]);
        let mut t1 = state[2].clone();
        t1.add_assign(&state[3]);

        let mut t2 = state[1].clone();
        t2.double();
        t2.add_assign(&t1);

        let mut t3 = state[3].clone();
        t3.double();
        t3.add_assign(&t0);

        let mut t4 = t1;
        t4.double();
        t4.double();
        t4.add_assign(&t3);

        let mut t5 = t0;
        t5.double();
        t5.double();
        t5.add_assign(&t2);

        let mut t6 = t3;
        t6.add_assign(&t5);

        let mut t7 = t2;
        t7.add_assign(&t4);

        state[0] = t6;
        state[1] = t5;
        state[2] = t7;
        state[3] = t4;
    }

    fn add_rc_in_place(&self, state: &mut [F], round: usize) {
        let rc = &self.params.round_constants[round];
        for (x, c) in state.iter_mut().zip(rc.iter()) {
            x.add_assign(c);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::babybear::BabyBear;
    use crate::fields::koalabear::KoalaBear;
    use crate::fields::mersenne31::Mersenne31;
    use crate::griffin::instances::*;

    fn check_small_prime_kat<F: FieldElement>(
        params: &Arc<GriffinParams<F>>,
        expected: &[u64],
        rounds: usize,
    ) {
        let permutation = Griffin::new(params);
        let input = (1..=permutation.get_t())
            .map(|value| F::from_u64(value as u64))
            .collect::<Vec<_>>();
        let expected = expected
            .iter()
            .copied()
            .map(F::from_u64)
            .collect::<Vec<_>>();
        assert_eq!(permutation.permutation(&input), expected);
        assert_eq!(params.rounds, rounds);
    }

    #[test]
    fn generated_small_prime_instances_match_reference_kats() {
        check_small_prime_kat::<BabyBear>(
            &GRIFFIN_BABYBEAR_16_PARAMS,
            &[
                1797865528, 1225257946, 790313438, 139513769, 1975768416, 1962642100, 546307734,
                1237972297, 1478591913, 980513569, 1733496131, 916946445, 273212616, 909656049,
                1427501869, 502580295,
            ],
            15,
        );
        check_small_prime_kat::<BabyBear>(
            &GRIFFIN_BABYBEAR_24_PARAMS,
            &[
                947136114, 1154062, 1375347716, 1532147275, 1897054674, 1523101362, 1052556190,
                644880643, 501346622, 903458917, 1940646043, 428971352, 835395889, 1911171217,
                942069770, 671069317, 1471865100, 1785099512, 1759805785, 1174582381, 1401983777,
                370744117, 599523798, 440971523,
            ],
            15,
        );
        check_small_prime_kat::<KoalaBear>(
            &GRIFFIN_KOALABEAR_16_PARAMS,
            &[
                579303170, 1946082004, 1207882601, 626619282, 674247378, 1187972934, 744390244,
                2085751129, 1724042636, 1756244318, 2068207835, 929706339, 875289086, 1518284839,
                1971355310, 1857739395,
            ],
            14,
        );
        check_small_prime_kat::<KoalaBear>(
            &GRIFFIN_KOALABEAR_24_PARAMS,
            &[
                1790550079, 1870580895, 1869198660, 242617202, 1249492959, 391173774, 772392565,
                318013841, 987231980, 1683296758, 293228204, 1385393336, 317417275, 1155311415,
                1782310590, 2027951434, 146970653, 562712538, 1698385074, 424004726, 676426835,
                1093692456, 1260937829, 1064636120,
            ],
            14,
        );
        check_small_prime_kat::<Mersenne31>(
            &GRIFFIN_MERSENNE31_16_PARAMS,
            &[
                23579611, 1195115056, 781712901, 1145417403, 1407753185, 26632645, 825054599,
                19569071, 210603208, 1606046957, 1924246711, 585422741, 302966780, 1098746637,
                38806779, 1806122786,
            ],
            15,
        );
        check_small_prime_kat::<Mersenne31>(
            &GRIFFIN_MERSENNE31_24_PARAMS,
            &[
                382894519, 1241143285, 753697801, 1681061897, 302492270, 988184361, 1924690703,
                2134670966, 61306371, 2139338212, 202392292, 1471733421, 730828355, 891447803,
                1336518768, 1402965022, 843463982, 282044900, 1232644666, 2031641104, 464021205,
                968657599, 334869526, 1524380121,
            ],
            15,
        );
    }

    fn assert_same_non_round_profile<F: FieldElement>(
        current: &GriffinParams<F>,
        original: &GriffinParams<F>,
    ) {
        assert_eq!(current.t, original.t);
        assert_eq!(current.d, original.d);
        assert_eq!(current.d_inv, original.d_inv);
        assert_eq!(
            &current.round_constants[..original.rounds - 1],
            original.round_constants
        );
    }

    #[test]
    fn paired_profiles_share_all_non_round_choices() {
        assert_same_non_round_profile(
            &GRIFFIN_BN254_3_CURRENT_PARAMS,
            &GRIFFIN_BN254_3_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &GRIFFIN_BLS12_381_3_CURRENT_PARAMS,
            &GRIFFIN_BLS12_381_3_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &GRIFFIN_BN254_4_CURRENT_PARAMS,
            &GRIFFIN_BN254_4_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &GRIFFIN_BLS12_381_4_CURRENT_PARAMS,
            &GRIFFIN_BLS12_381_4_ORIGINAL_PARAMS,
        );
    }
}
