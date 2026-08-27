use crate::fields::FieldElement;
use std::sync::Arc;

/// GMiMCHash2 — GMiMC^π_erf2 permutation as specified in §5.1.2 of the SoK paper.
///
/// Differences from the original GMiMC_erf:
///   1. S-box exponent α is a power of two (the legacy instances use 2; the
///      recommended large-field t=4 profile uses 8).
///   2. An Input/Output matrix `M_IO` is applied before the first round and after the
///      last round. For 256-bit fields (c=1) it is the identity; for smaller fields it is
///      a sparse circulant matrix (see §5.1.2, page 25).
///   3. Round numbers are taken from Table 9 of the SoK paper.
///
/// Permutation structure:
///   π(x) = M_IO ∘ R^(R) ∘ … ∘ R^(1) (M_IO · x)
///
/// Round function R^(i) (Figure 5):
///   state[0]  += rc[i]                 // constant addition
///   pow        = state[0]²             // square S-box
///   state[1..] += pow                  // distribute
///   rotate_right(1)  (all rounds except the last)

#[derive(Clone, Debug)]
pub struct Gmimc2Params<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) alpha: u64,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<F>,
}

impl<F: FieldElement> Gmimc2Params<F> {
    pub fn new(t: usize, rounds: usize, round_constants: &[F]) -> Self {
        Self::new_with_alpha(t, 2, rounds, round_constants)
    }

    pub fn new_with_alpha(t: usize, alpha: u64, rounds: usize, round_constants: &[F]) -> Self {
        assert!(alpha >= 2 && alpha.is_power_of_two());
        assert_eq!(round_constants.len(), rounds);
        Gmimc2Params {
            t,
            alpha,
            rounds,
            round_constants: round_constants.to_owned(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Gmimc2<F: FieldElement> {
    pub(crate) params: Arc<Gmimc2Params<F>>,
}

impl<F: FieldElement> Gmimc2<F> {
    pub fn new(params: &Arc<Gmimc2Params<F>>) -> Self {
        Gmimc2 {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    // ── Input/Output matrix (§5.1.2, page 25) ──
    //
    //  c=1 (256-bit)           → identity
    //  t∈{8,16}, c=r=t/2      → circ(1,0,…,0,2,0,…,0)   with 2 at offset t/2
    //  t∈{12,24}, c=t/3       → circ(1,0,…,0,2,0,…,0,2,0,…,0)  with 2s at t/3 and t/2
    //  t=4                    → identity (implicit from c=1 rule)

    fn apply_circulant(state: &[F]) -> Vec<F> {
        let t = state.len();
        let offsets: &[usize] = match t {
            // 256-bit fields: identity
            2 | 3 | 4 | 5 => return state.to_vec(),
            // t∈{8,16}: 2 at t/2
            8 => &[4],
            16 => &[8],
            // t∈{12,24}: 2 at t/3 and t/2
            12 => &[4, 6],
            24 => &[8, 12],
            // fallback
            _ => &[t / 2],
        };
        let mut result = state.to_vec();
        for (i, ri) in result.iter_mut().enumerate() {
            for &off in offsets {
                let mut term = state[(i + off) % t].clone();
                term.double();
                ri.add_assign(&term);
            }
        }
        result
    }

    // ── S-box: (x + rc)^alpha ──

    #[inline(always)]
    fn sbox_p(&self, x: &F) -> F {
        x.pow_u64(self.params.alpha)
    }

    // ── Permutation ──

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        let r = self.params.rounds;
        assert_eq!(input.len(), t);
        if r == 0 {
            return input.to_vec();
        }

        let mut state = Self::apply_circulant(input);

        // The GMiMC2 constant is part of branch zero and travels with that
        // branch after the cyclic left shift; it is not merely an S-box input.
        for round in 0..r {
            state[0].add_assign(&self.params.round_constants[round]);
            let pow = self.sbox_p(&state[0]);
            for el in state.iter_mut().skip(1) {
                el.add_assign(&pow);
            }
            state.rotate_left(1);
        }
        Self::apply_circulant(&state)
    }
}

#[cfg(test)]
mod tests {
    use super::{Gmimc2, Gmimc2Params};
    use crate::fields::babybear::BabyBear;
    use crate::fields::bls12_381::Bls12_381;
    use crate::fields::bn254::Bn254;
    use crate::fields::goldilocks::Goldilocks;
    use crate::fields::koalabear::KoalaBear;
    use crate::fields::mersenne31::Mersenne31;
    use crate::fields::FieldElement;
    use crate::gmimc2::instances::*;
    use std::sync::Arc;

    fn fingerprint<F: FieldElement>(values: &[F]) -> F {
        let multiplier = F::from_u64(257);
        values.iter().fold(F::zero(), |mut acc, value| {
            acc.mul_assign(&multiplier);
            acc.add_assign(value);
            acc
        })
    }

    fn check_profile<F: FieldElement>(
        params: &Arc<Gmimc2Params<F>>,
        t: usize,
        alpha: u64,
        rounds: usize,
        expected_fingerprint: F,
    ) {
        assert_eq!(params.t, t);
        assert_eq!(params.alpha, alpha);
        assert_eq!(params.rounds, rounds);
        assert_eq!(params.round_constants.len(), rounds);
        let input = (1..=t)
            .map(|value| F::from_u64(value as u64))
            .collect::<Vec<_>>();
        let output = Gmimc2::new(params).permutation(&input);
        assert_eq!(fingerprint(&output), expected_fingerprint);
    }

    #[test]
    fn source_selected_profiles_and_kats_are_exact() {
        check_profile(
            &GMIMC2_BLS12_381_3_PARAMS,
            3,
            8,
            63,
            Bls12_381::from_hex("32647337cdb393bbef50ac22d0b4a3ca617f63c8ff3979f51f86b0c7a2e2fa0b")
                .unwrap(),
        );
        check_profile(
            &GMIMC2_BN254_3_PARAMS,
            3,
            8,
            63,
            Bn254::from_hex("20ab033578ed69911a77bd9db69baa9f6de0b7f2bdb43abcb847d33350876eef")
                .unwrap(),
        );
        check_profile(
            &GMIMC2_GOLDILOCKS_8_PARAMS,
            8,
            4,
            88,
            Goldilocks::from_u64(17_166_913_452_726_183_267),
        );
        check_profile(
            &GMIMC2_GOLDILOCKS_12_PARAMS,
            12,
            4,
            96,
            Goldilocks::from_u64(4_426_234_654_506_097_107),
        );
        check_profile(
            &GMIMC2_BABYBEAR_16_PARAMS,
            16,
            2,
            176,
            BabyBear::from_u64(659_549_587),
        );
        check_profile(
            &GMIMC2_BABYBEAR_24_PARAMS,
            24,
            2,
            264,
            BabyBear::from_u64(1_243_511_552),
        );
        check_profile(
            &GMIMC2_KOALABEAR_16_PARAMS,
            16,
            2,
            176,
            KoalaBear::from_u64(834_174_982),
        );
        check_profile(
            &GMIMC2_KOALABEAR_24_PARAMS,
            24,
            2,
            264,
            KoalaBear::from_u64(864_370_765),
        );
        check_profile(
            &GMIMC2_MERSENNE31_16_PARAMS,
            16,
            2,
            176,
            Mersenne31::from_u64(1_666_050_320),
        );
        check_profile(
            &GMIMC2_MERSENNE31_24_PARAMS,
            24,
            2,
            264,
            Mersenne31::from_u64(1_870_415_535),
        );
    }
}
