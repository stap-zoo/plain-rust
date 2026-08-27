use crate::fields::FieldElement;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct RescuePrimeParams<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) d: u64,
    pub(crate) d_inv: [u64; 4],
    pub(crate) rounds: usize,
    pub(crate) mds: Vec<Vec<F>>,
    pub(crate) round_constants: Vec<Vec<F>>, // [round_idx][state_idx], round_idx in [0, 2 * rounds)
}

impl<F: FieldElement> RescuePrimeParams<F> {
    pub fn new(
        t: usize,
        d: u64,
        d_inv_words: [u64; 4],
        rounds: usize,
        mds: &[Vec<F>],
        round_constants: &[Vec<F>],
    ) -> Self {
        assert!(d == 3 || d == 5 || d == 7);
        assert_eq!(mds.len(), t);
        for row in mds {
            assert_eq!(row.len(), t);
        }
        assert_eq!(round_constants.len(), 2 * rounds);
        for rc in round_constants {
            assert_eq!(rc.len(), t);
        }

        RescuePrimeParams {
            t,
            d,
            d_inv: d_inv_words,
            rounds,
            mds: mds.to_owned(),
            round_constants: round_constants.to_owned(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct RescuePrime<F: FieldElement> {
    pub(crate) params: Arc<RescuePrimeParams<F>>,
}

impl<F: FieldElement> RescuePrime<F> {
    pub fn new(params: &Arc<RescuePrimeParams<F>>) -> Self {
        RescuePrime {
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
        for r in 0..self.params.rounds {
            for x in state.iter_mut() {
                *x = self.sbox_p(x);
            }
            self.affine_round(&mut state, 2 * r);

            for x in state.iter_mut() {
                *x = self.sbox_p_inv(x);
            }
            self.affine_round(&mut state, 2 * r + 1);
        }
        state
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

    fn sbox_p_inv(&self, input: &F) -> F {
        input.pow_words_le(&self.params.d_inv)
    }

    fn affine_round(&self, state: &mut [F], round: usize) {
        self.matmul_in_place(state, &self.params.mds);
        self.add_rc_in_place(state, round);
    }

    fn matmul_in_place(&self, state: &mut [F], mat: &[Vec<F>]) {
        let t = state.len();
        let mut out = vec![F::zero(); t];
        for row in 0..t {
            for (col, val) in state.iter().enumerate().take(t) {
                let mut tmp = mat[row][col].clone();
                tmp.mul_assign(val);
                out[row].add_assign(&tmp);
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::babybear::BabyBear;
    use crate::fields::koalabear::KoalaBear;
    use crate::fields::mersenne31::Mersenne31;
    use crate::rescueprime::instances::*;

    fn check_small_prime_kat<F: FieldElement>(
        params: &Arc<RescuePrimeParams<F>>,
        expected: &[u64],
    ) {
        let permutation = RescuePrime::new(params);
        let input = (1..=permutation.get_t())
            .map(|value| F::from_u64(value as u64))
            .collect::<Vec<_>>();
        let expected = expected
            .iter()
            .copied()
            .map(F::from_u64)
            .collect::<Vec<_>>();
        assert_eq!(permutation.permutation(&input), expected);
        assert_eq!(params.rounds, 8);
    }

    #[test]
    fn generated_small_prime_instances_match_reference_kats() {
        check_small_prime_kat::<BabyBear>(
            &RESCUE_PRIME_BABYBEAR_16_PARAMS,
            &[
                148520378, 750979325, 1742931594, 938695099, 1876875392, 61759309, 775195108,
                1520242022, 561754184, 1599021303, 836884156, 1252837894, 1547265626, 1418428354,
                1235082754, 779920601,
            ],
        );
        check_small_prime_kat::<BabyBear>(
            &RESCUE_PRIME_BABYBEAR_24_PARAMS,
            &[
                9377853, 2010981566, 1135267071, 867286365, 701292781, 1541004182, 1871056350,
                1182766233, 1591900104, 1352961775, 299933196, 726875324, 528894994, 1788642950,
                1771855496, 1176234109, 1864024610, 593734090, 1682519838, 1695591087, 226317897,
                1947486960, 491944121, 339704787,
            ],
        );
        check_small_prime_kat::<KoalaBear>(
            &RESCUE_PRIME_KOALABEAR_16_PARAMS,
            &[
                406688936, 2098855104, 862866808, 71515720, 1141886882, 1095014288, 587056277,
                1537902750, 1224812489, 1617539216, 496543356, 230424620, 1109727730, 870182440,
                1538048363, 519566810,
            ],
        );
        check_small_prime_kat::<KoalaBear>(
            &RESCUE_PRIME_KOALABEAR_24_PARAMS,
            &[
                1515762978, 1772160492, 1939596788, 1295816424, 45118571, 502770126, 1095054859,
                459091376, 1391031976, 2061166925, 952270003, 892175529, 678766284, 1906994494,
                1793462039, 1412627102, 1650261071, 1875798583, 722494996, 383285076, 426379969,
                1852972083, 1749085800, 1780070173,
            ],
        );
        check_small_prime_kat::<Mersenne31>(
            &RESCUE_PRIME_MERSENNE31_16_PARAMS,
            &[
                1136589846, 1908178910, 1312842829, 1117634965, 640577991, 325267093, 930267757,
                1207850025, 568429631, 1169853887, 1766710223, 796615643, 575313684, 1171985373,
                530691354, 948914467,
            ],
        );
        check_small_prime_kat::<Mersenne31>(
            &RESCUE_PRIME_MERSENNE31_24_PARAMS,
            &[
                80405035, 374334217, 1589360205, 2097587048, 976354558, 2101374562, 611332052,
                1564257281, 1450168624, 1261350120, 94483386, 172224239, 855242256, 553492948,
                1671177689, 793656732, 82109724, 1745335452, 454738199, 208291330, 1408771080,
                266949708, 831288684, 222181969,
            ],
        );
    }

    fn assert_same_non_round_profile<F: FieldElement>(
        current: &RescuePrimeParams<F>,
        original: &RescuePrimeParams<F>,
    ) {
        assert_eq!(current.t, original.t);
        assert_eq!(current.d, original.d);
        assert_eq!(current.d_inv, original.d_inv);
        assert_eq!(current.mds, original.mds);
        assert_eq!(
            &current.round_constants[..2 * original.rounds],
            original.round_constants
        );
    }

    #[test]
    fn paired_profiles_share_all_non_round_choices() {
        assert_same_non_round_profile(
            &RESCUE_PRIME_BN254_3_CURRENT_PARAMS,
            &RESCUE_PRIME_BN254_3_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &RESCUE_PRIME_BLS12_381_3_CURRENT_PARAMS,
            &RESCUE_PRIME_BLS12_381_3_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &RESCUE_PRIME_BN254_4_CURRENT_PARAMS,
            &RESCUE_PRIME_BN254_4_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &RESCUE_PRIME_BLS12_381_4_CURRENT_PARAMS,
            &RESCUE_PRIME_BLS12_381_4_ORIGINAL_PARAMS,
        );
    }
}
