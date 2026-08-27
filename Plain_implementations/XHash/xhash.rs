//Plain implementation of the four formal DCC 2026 XHash profiles, rather than the earlier three-round draft retained elsewhere in this repository.
//XHash8/12 use a 12-word Goldilocks state and XHash16/24 use a 24-word Mersenne31 state; the instance name counts active inverse S-box lanes, not state width.
//Every profile uses six rounds and constants derived with the authors' official SHAKE256 procedure, with no extra final linear layer or sponge-level processing.

use crate::fields::FieldElement;
use std::sync::Arc;

/// Number of rounds in the formal DCC 2026 XHash instances.
pub const XHASH_ROUNDS: usize = 6;

/// Field-specific XHash profile.
///
/// Besides selecting the base-field exponent, a profile fixes the cubic
/// extension used by the extension-field S-box:
///
/// - Goldilocks: `F_p[u] / (u^3 - u - 1)` and `alpha = 7`;
/// - Mersenne31: `F_p[u] / (u^3 + 2)` and `alpha = 5`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum XHashProfile {
    Goldilocks,
    Mersenne31,
}

/// Choice of the inverse-power layer in an XHash round.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InverseLayer {
    /// Apply the inverse power only at lanes `j` for which `j % 3 != 1`.
    Partial,
    /// Apply the inverse power to every lane.
    Full,
}

impl InverseLayer {
    #[inline(always)]
    fn is_active(self, lane: usize) -> bool {
        self == Self::Full || lane % 3 != 1
    }
}

/// Parameters for one of the four formal XHash instances.
#[derive(Clone, Debug)]
pub struct XHashParams<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) capacity: usize,
    pub(crate) rounds: usize,
    pub(crate) alpha: u64,
    pub(crate) alpha_inv: u64,
    pub(crate) profile: XHashProfile,
    pub(crate) inverse_layer: InverseLayer,
    pub(crate) round_constants: Vec<F>,
    pub(crate) mds_matrix: Vec<Vec<F>>,
}

impl<F: FieldElement> XHashParams<F> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t: usize,
        capacity: usize,
        rounds: usize,
        alpha: u64,
        alpha_inv: u64,
        profile: XHashProfile,
        inverse_layer: InverseLayer,
        round_constants: &[F],
        mds_matrix: Vec<Vec<F>>,
    ) -> Self {
        assert!(matches!(t, 12 | 24), "XHash state width must be 12 or 24");
        assert!(
            capacity < t,
            "XHash capacity must be smaller than its state"
        );
        assert_eq!(
            rounds, XHASH_ROUNDS,
            "formal XHash instances use exactly six rounds"
        );
        assert_eq!(
            round_constants.len(),
            3 * t * rounds,
            "XHash needs three width-sized constant blocks per round"
        );
        assert_eq!(mds_matrix.len(), t, "XHash MDS matrix has wrong height");
        assert!(
            mds_matrix.iter().all(|row| row.len() == t),
            "XHash MDS matrix must be square"
        );

        match profile {
            XHashProfile::Goldilocks => {
                assert_eq!(t, 12, "Goldilocks XHash uses a 12-word state");
                assert_eq!(capacity, 4, "Goldilocks XHash uses capacity four");
                assert_eq!(alpha, 7, "Goldilocks XHash uses alpha = 7");
                assert_eq!(
                    alpha_inv, 10_540_996_611_094_048_183,
                    "wrong inverse exponent for Goldilocks alpha = 7"
                );
            }
            XHashProfile::Mersenne31 => {
                assert_eq!(t, 24, "M31 XHash uses a 24-word state");
                assert_eq!(capacity, 8, "M31 XHash uses capacity eight");
                assert_eq!(alpha, 5, "M31 XHash uses alpha = 5");
                assert_eq!(
                    alpha_inv, 1_717_986_917,
                    "wrong inverse exponent for M31 alpha = 5"
                );
            }
        }

        Self {
            t,
            capacity,
            rounds,
            alpha,
            alpha_inv,
            profile,
            inverse_layer,
            round_constants: round_constants.to_vec(),
            mds_matrix,
        }
    }

    pub fn get_t(&self) -> usize {
        self.t
    }

    pub fn get_rate(&self) -> usize {
        self.t - self.capacity
    }

    pub fn get_capacity(&self) -> usize {
        self.capacity
    }

    pub fn get_rounds(&self) -> usize {
        self.rounds
    }

    pub fn profile(&self) -> XHashProfile {
        self.profile
    }

    pub fn inverse_layer(&self) -> InverseLayer {
        self.inverse_layer
    }

    pub fn active_inverse_sboxes(&self) -> usize {
        match self.inverse_layer {
            InverseLayer::Partial => (0..self.t).filter(|&j| j % 3 != 1).count(),
            InverseLayer::Full => self.t,
        }
    }

    pub fn round_constants(&self) -> &[F] {
        &self.round_constants
    }

    pub fn mds_matrix(&self) -> &[Vec<F>] {
        &self.mds_matrix
    }
}

/// The formal XHash permutation.
#[derive(Clone, Debug)]
pub struct XHash<F: FieldElement> {
    params: Arc<XHashParams<F>>,
}

impl<F: FieldElement> XHash<F> {
    pub fn new(params: &Arc<XHashParams<F>>) -> Self {
        Self {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn get_rounds(&self) -> usize {
        self.params.rounds
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let mut state = input.to_vec();
        self.permutation_in_place(&mut state);
        state
    }

    /// Evaluate the DCC 2026 round function in place.
    ///
    /// Each round is
    /// `C_F -> MDS -> x^alpha -> MDS -> C_B -> x^(1/alpha)
    ///  -> C_E -> (Fp3 x^alpha)`.
    /// There are six rounds and no additional final linear step.
    pub fn permutation_in_place(&self, state: &mut [F]) {
        self.validate_state(state);

        for round in 0..self.params.rounds {
            self.add_constants(state, round, 0);
            self.mds_layer_in_place(state);
            self.base_power_layer(state);

            self.mds_layer_in_place(state);
            self.add_constants(state, round, 1);
            self.inverse_power_layer(state);

            self.add_constants(state, round, 2);
            self.extension_power_layer(state);
        }
    }

    /// Alias matching the in-place naming used by permutation libraries.
    pub fn permute_in_place(&self, state: &mut [F]) {
        self.permutation_in_place(state);
    }

    /// Literal allocating implementation retained as an audit oracle.
    ///
    /// This follows the written round definition directly and uses generic
    /// square-and-multiply in the cubic field. The in-place path uses the
    /// fixed `x^5`/`x^7` addition chains instead.
    pub fn permutation_reference(&self, input: &[F]) -> Vec<F> {
        self.validate_state(input);
        let mut state = input.to_vec();

        for round in 0..self.params.rounds {
            add_constant_block(&mut state, self.constant_block(round, 0));
            state = mds_matmul(&state, &self.params.mds_matrix);
            for value in &mut state {
                *value = value.pow_u64(self.params.alpha);
            }

            state = mds_matmul(&state, &self.params.mds_matrix);
            add_constant_block(&mut state, self.constant_block(round, 1));
            for (lane, value) in state.iter_mut().enumerate() {
                if self.params.inverse_layer.is_active(lane) {
                    *value = value.pow_u64(self.params.alpha_inv);
                }
            }

            add_constant_block(&mut state, self.constant_block(round, 2));
            for chunk in state.chunks_exact_mut(3) {
                let value = [chunk[0].clone(), chunk[1].clone(), chunk[2].clone()];
                let raised = cubic_pow(value, self.params.alpha, self.params.profile);
                chunk.clone_from_slice(&raised);
            }
        }

        state
    }

    fn validate_state(&self, state: &[F]) {
        assert_eq!(
            state.len(),
            self.params.t,
            "XHash input width does not match its parameters"
        );
    }

    fn constant_block(&self, round: usize, step: usize) -> &[F] {
        let start = (3 * round + step) * self.params.t;
        &self.params.round_constants[start..start + self.params.t]
    }

    #[inline(always)]
    fn add_constants(&self, state: &mut [F], round: usize, step: usize) {
        add_constant_block(state, self.constant_block(round, step));
    }

    #[inline(always)]
    fn mds_layer_in_place(&self, state: &mut [F]) {
        // Both formal widths fit in this fixed scratch buffer. Keeping it on
        // the stack avoids two heap allocations per XHash round.
        let mut output: [F; 24] = core::array::from_fn(|_| F::zero());
        for (row, output_value) in output.iter_mut().enumerate().take(state.len()) {
            for (coefficient, input_value) in self.params.mds_matrix[row].iter().zip(state.iter()) {
                let mut product = coefficient.clone();
                product.mul_assign(input_value);
                output_value.add_assign(&product);
            }
        }
        state.clone_from_slice(&output[..state.len()]);
    }

    #[inline(always)]
    fn base_power_layer(&self, state: &mut [F]) {
        for value in state {
            let original = value.clone();
            let mut square = original.clone();
            square.square();
            let mut fourth = square.clone();
            fourth.square();

            match self.params.alpha {
                5 => {
                    fourth.mul_assign(&original);
                    *value = fourth;
                }
                7 => {
                    fourth.mul_assign(&square);
                    fourth.mul_assign(&original);
                    *value = fourth;
                }
                _ => unreachable!("formal XHash profiles use alpha 5 or 7"),
            }
        }
    }

    #[inline(always)]
    fn inverse_power_layer(&self, state: &mut [F]) {
        for (lane, value) in state.iter_mut().enumerate() {
            if self.params.inverse_layer.is_active(lane) {
                *value = value.pow_u64(self.params.alpha_inv);
            }
        }
    }

    #[inline(always)]
    fn extension_power_layer(&self, state: &mut [F]) {
        for chunk in state.chunks_exact_mut(3) {
            let value = [chunk[0].clone(), chunk[1].clone(), chunk[2].clone()];
            let raised = match self.params.alpha {
                5 => cubic_pow5(value, self.params.profile),
                7 => cubic_pow7(value, self.params.profile),
                _ => unreachable!("formal XHash profiles use alpha 5 or 7"),
            };
            chunk.clone_from_slice(&raised);
        }
    }
}

#[inline(always)]
fn add_constant_block<F: FieldElement>(state: &mut [F], constants: &[F]) {
    debug_assert_eq!(state.len(), constants.len());
    for (value, constant) in state.iter_mut().zip(constants) {
        value.add_assign(constant);
    }
}

fn mds_matmul<F: FieldElement>(state: &[F], matrix: &[Vec<F>]) -> Vec<F> {
    let mut output = vec![F::zero(); state.len()];
    for (row, output_value) in output.iter_mut().enumerate() {
        for (coefficient, input_value) in matrix[row].iter().zip(state) {
            let mut product = coefficient.clone();
            product.mul_assign(input_value);
            output_value.add_assign(&product);
        }
    }
    output
}

#[inline(always)]
fn cubic_pow5<F: FieldElement>(value: [F; 3], profile: XHashProfile) -> [F; 3] {
    let square = cubic_mul(&value, &value, profile);
    let fourth = cubic_mul(&square, &square, profile);
    cubic_mul(&fourth, &value, profile)
}

#[inline(always)]
fn cubic_pow7<F: FieldElement>(value: [F; 3], profile: XHashProfile) -> [F; 3] {
    let square = cubic_mul(&value, &value, profile);
    let fourth = cubic_mul(&square, &square, profile);
    let sixth = cubic_mul(&fourth, &square, profile);
    cubic_mul(&sixth, &value, profile)
}

fn cubic_pow<F: FieldElement>(
    mut base: [F; 3],
    mut exponent: u64,
    profile: XHashProfile,
) -> [F; 3] {
    let mut result = [F::one(), F::zero(), F::zero()];
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = cubic_mul(&result, &base, profile);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = cubic_mul(&base, &base, profile);
        }
    }
    result
}

/// Multiply in the profile's cubic extension, represented in the basis
/// `(1, u, u^2)`.
fn cubic_mul<F: FieldElement>(left: &[F; 3], right: &[F; 3], profile: XHashProfile) -> [F; 3] {
    let d0 = product(&left[0], &right[0]);

    let mut d1 = product(&left[0], &right[1]);
    d1.add_assign(&product(&left[1], &right[0]));

    let mut d2 = product(&left[0], &right[2]);
    d2.add_assign(&product(&left[1], &right[1]));
    d2.add_assign(&product(&left[2], &right[0]));

    let mut d3 = product(&left[1], &right[2]);
    d3.add_assign(&product(&left[2], &right[1]));

    let d4 = product(&left[2], &right[2]);

    match profile {
        // u^3 = u + 1 and u^4 = u^2 + u.
        XHashProfile::Goldilocks => {
            let mut c0 = d0;
            c0.add_assign(&d3);

            let mut c1 = d1;
            c1.add_assign(&d3);
            c1.add_assign(&d4);

            let mut c2 = d2;
            c2.add_assign(&d4);
            [c0, c1, c2]
        }
        // u^3 = -2 and u^4 = -2u.
        XHashProfile::Mersenne31 => {
            let mut twice_d3 = d3;
            twice_d3.double();
            let mut c0 = d0;
            c0.sub_assign(&twice_d3);

            let mut twice_d4 = d4;
            twice_d4.double();
            let mut c1 = d1;
            c1.sub_assign(&twice_d4);

            [c0, c1, d2]
        }
    }
}

#[inline(always)]
fn product<F: FieldElement>(left: &F, right: &F) -> F {
    let mut result = left.clone();
    result.mul_assign(right);
    result
}

#[cfg(test)]
mod tests {
    use super::super::instances::{
        XHASH12_GOLDILOCKS_PARAMS, XHASH16_M31_PARAMS, XHASH24_M31_PARAMS, XHASH8_GOLDILOCKS_PARAMS,
    };
    use super::{cubic_mul, InverseLayer, XHash, XHashParams, XHashProfile, XHASH_ROUNDS};
    use crate::fields::goldilocks::Goldilocks;
    use crate::fields::mersenne31::Mersenne31;
    use crate::fields::FieldElement;

    fn deterministic_goldilocks_input(seed: u64) -> Vec<Goldilocks> {
        let mut value = seed;
        (0..12)
            .map(|_| {
                value = value
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                Goldilocks::from_u64(value)
            })
            .collect()
    }

    fn deterministic_m31_input(seed: u32) -> Vec<Mersenne31> {
        let mut value = seed;
        (0..24)
            .map(|_| {
                value = value.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                Mersenne31::from_u32(value % Mersenne31::MODULUS)
            })
            .collect()
    }

    #[test]
    fn optimized_matches_literal_reference_for_all_formal_instances() {
        let gold8 = XHash::new(&XHASH8_GOLDILOCKS_PARAMS);
        let gold12 = XHash::new(&XHASH12_GOLDILOCKS_PARAMS);
        let m31_16 = XHash::new(&XHASH16_M31_PARAMS);
        let m31_24 = XHash::new(&XHASH24_M31_PARAMS);

        for seed in [0, 1, 42, u64::MAX] {
            let input = deterministic_goldilocks_input(seed);
            assert_eq!(
                gold8.permutation(&input),
                gold8.permutation_reference(&input)
            );
            assert_eq!(
                gold12.permutation(&input),
                gold12.permutation_reference(&input)
            );
        }

        for seed in [0, 1, 42, u32::MAX] {
            let input = deterministic_m31_input(seed);
            assert_eq!(
                m31_16.permutation(&input),
                m31_16.permutation_reference(&input)
            );
            assert_eq!(
                m31_24.permutation(&input),
                m31_24.permutation_reference(&input)
            );
        }
    }

    #[test]
    fn formal_instance_shapes_and_partial_mask_are_exact() {
        assert_eq!(XHASH8_GOLDILOCKS_PARAMS.get_t(), 12);
        assert_eq!(XHASH8_GOLDILOCKS_PARAMS.rounds, 6);
        assert_eq!(XHASH8_GOLDILOCKS_PARAMS.alpha, 7);
        assert_eq!(XHASH8_GOLDILOCKS_PARAMS.get_rate(), 8);
        assert_eq!(XHASH8_GOLDILOCKS_PARAMS.active_inverse_sboxes(), 8);
        assert_eq!(
            XHASH8_GOLDILOCKS_PARAMS.inverse_layer(),
            InverseLayer::Partial
        );

        assert_eq!(XHASH12_GOLDILOCKS_PARAMS.get_t(), 12);
        assert_eq!(XHASH12_GOLDILOCKS_PARAMS.rounds, 6);
        assert_eq!(XHASH12_GOLDILOCKS_PARAMS.alpha, 7);
        assert_eq!(XHASH12_GOLDILOCKS_PARAMS.active_inverse_sboxes(), 12);
        assert_eq!(
            XHASH12_GOLDILOCKS_PARAMS.inverse_layer(),
            InverseLayer::Full
        );

        assert_eq!(XHASH16_M31_PARAMS.get_t(), 24);
        assert_eq!(XHASH16_M31_PARAMS.rounds, 6);
        assert_eq!(XHASH16_M31_PARAMS.alpha, 5);
        assert_eq!(XHASH16_M31_PARAMS.get_rate(), 16);
        assert_eq!(XHASH16_M31_PARAMS.active_inverse_sboxes(), 16);
        assert_eq!(XHASH16_M31_PARAMS.inverse_layer(), InverseLayer::Partial);

        assert_eq!(XHASH24_M31_PARAMS.get_t(), 24);
        assert_eq!(XHASH24_M31_PARAMS.rounds, 6);
        assert_eq!(XHASH24_M31_PARAMS.alpha, 5);
        assert_eq!(XHASH24_M31_PARAMS.active_inverse_sboxes(), 24);
        assert_eq!(XHASH24_M31_PARAMS.inverse_layer(), InverseLayer::Full);
    }

    #[test]
    fn inverse_exponents_round_trip_the_forward_sboxes() {
        for value in [0, 1, 2, 7, Goldilocks::MODULUS - 1] {
            let x = Goldilocks::from_u64(value);
            assert_eq!(x.pow_u64(7).pow_u64(10_540_996_611_094_048_183), x);
        }
        for value in [0, 1, 2, 7, Mersenne31::MODULUS - 1] {
            let x = Mersenne31::from_u32(value);
            assert_eq!(x.pow_u64(5).pow_u64(1_717_986_917), x);
        }
    }

    #[test]
    fn cubic_reductions_match_the_formal_polynomials() {
        let zero_g = Goldilocks::zero();
        let one_g = Goldilocks::one();
        let u2_g = [zero_g, zero_g, one_g];
        let u_g = [zero_g, one_g, zero_g];
        assert_eq!(
            cubic_mul(&u2_g, &u_g, XHashProfile::Goldilocks),
            [one_g, one_g, zero_g]
        );

        let zero_m = Mersenne31::zero();
        let one_m = Mersenne31::one();
        let u2_m = [zero_m, zero_m, one_m];
        let u_m = [zero_m, one_m, zero_m];
        let minus_two = Mersenne31::from_u32(Mersenne31::MODULUS - 2);
        assert_eq!(
            cubic_mul(&u2_m, &u_m, XHashProfile::Mersenne31),
            [minus_two, zero_m, zero_m]
        );
    }

    /// Regression vectors for the six-round DCC 2026 profiles. They are not
    /// author-published KATs.
    #[test]
    fn independent_literal_reference_regression_vectors() {
        let gold_input: Vec<_> = (0..12).map(Goldilocks::from_u64).collect();
        let gold8_expected = [
            1_634_901_000_479_886_711,
            15_164_765_201_519_446_800,
            9_857_823_976_508_228_082,
            8_614_586_736_171_941_436,
            2_333_463_272_575_837_858,
            3_753_477_033_508_253_468,
            1_374_120_117_442_903_979,
            10_449_210_459_676_559_474,
            15_721_722_583_705_381_769,
            8_397_393_178_382_693_018,
            4_301_232_263_769_248_661,
            12_061_294_292_042_938_491,
        ];
        let gold12_expected = [
            3_376_642_126_165_865_447,
            3_246_235_210_200_447_438,
            7_368_853_067_160_481_556,
            5_002_761_933_388_891_071,
            9_609_111_108_636_321_968,
            17_758_540_707_947_177_461,
            16_337_707_017_830_674_577,
            2_257_955_953_815_099_177,
            5_232_398_019_987_582_873,
            12_164_611_353_555_196_203,
            5_832_817_359_388_281_458,
            15_900_478_155_810_051_674,
        ];
        let gold8_output = XHash::new(&XHASH8_GOLDILOCKS_PARAMS).permutation(&gold_input);
        let gold12_output = XHash::new(&XHASH12_GOLDILOCKS_PARAMS).permutation(&gold_input);
        assert_eq!(
            gold8_output
                .iter()
                .map(Goldilocks::to_u64)
                .collect::<Vec<_>>(),
            gold8_expected
        );
        assert_eq!(
            gold12_output
                .iter()
                .map(Goldilocks::to_u64)
                .collect::<Vec<_>>(),
            gold12_expected
        );

        let m31_input: Vec<_> = (0..24).map(Mersenne31::from_u32).collect();
        let m31_16_expected = [
            1_960_347_588,
            2_004_748_893,
            1_804_431_369,
            680_660_975,
            478_119_374,
            336_314_932,
            1_910_855_626,
            1_005_508_631,
            173_940_359,
            161_851_462,
            737_454_016,
            1_163_478_887,
            402_641_555,
            848_090_616,
            1_977_812_914,
            2_083_293_475,
            1_105_783_266,
            1_564_923_427,
            1_365_615_329,
            1_788_231_244,
            1_362_678_695,
            1_421_636_390,
            1_287_783_531,
            1_650_919_123,
        ];
        let m31_24_expected = [
            1_873_910_972,
            1_468_328_915,
            1_618_512_264,
            1_608_683_774,
            620_243_663,
            367_691_501,
            1_754_439_963,
            1_104_575_626,
            1_244_486_860,
            4_309_064,
            106_482_674,
            141_153_809,
            1_465_142_453,
            1_206_840_484,
            679_251_693,
            2_026_443_288,
            926_247_702,
            2_132_782_589,
            524_206_233,
            323_757_333,
            1_561_400_045,
            1_254_747_993,
            1_629_327_218,
            515_021_800,
        ];
        let m31_16_output = XHash::new(&XHASH16_M31_PARAMS).permutation(&m31_input);
        let m31_24_output = XHash::new(&XHASH24_M31_PARAMS).permutation(&m31_input);
        assert_eq!(
            m31_16_output
                .iter()
                .map(Mersenne31::to_u32)
                .collect::<Vec<_>>(),
            m31_16_expected
        );
        assert_eq!(
            m31_24_output
                .iter()
                .map(Mersenne31::to_u32)
                .collect::<Vec<_>>(),
            m31_24_expected
        );
    }

    #[test]
    #[should_panic(expected = "three width-sized constant blocks")]
    fn constructor_rejects_wrong_constant_count() {
        XHashParams::new(
            12,
            4,
            XHASH_ROUNDS,
            7,
            10_540_996_611_094_048_183,
            XHashProfile::Goldilocks,
            InverseLayer::Partial,
            &[],
            XHASH8_GOLDILOCKS_PARAMS.mds_matrix.clone(),
        );
    }

    #[test]
    #[should_panic(expected = "exactly six rounds")]
    fn constructor_rejects_wrong_round_count() {
        XHashParams::new(
            12,
            4,
            XHASH_ROUNDS - 1,
            7,
            10_540_996_611_094_048_183,
            XHashProfile::Goldilocks,
            InverseLayer::Partial,
            &XHASH8_GOLDILOCKS_PARAMS.round_constants,
            XHASH8_GOLDILOCKS_PARAMS.mds_matrix.clone(),
        );
    }
}
