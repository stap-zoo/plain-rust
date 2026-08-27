//Plain implementation of the formal untweaked pSquare-hash permutation over the 31-bit fields.
//It supports the paper's Mersenne31, KoalaBear, and BabyBear `t = 16/24`, capacity-8, 52-round instances, with constants and test vectors adapted from the Plonky3 implementations.
//The module provides both a literal reference path and an allocation-free double-round path; sponge padding and domain separation are outside its scope.

use crate::fields::FieldElement;
use std::sync::Arc;

/// Number of Feistel rounds in the formal untweaked pSquare-hash instances.
pub const PSQUAREHASH_ROUNDS: usize = 52;

const SUPPORTED_WIDTHS: [usize; 2] = [16, 24];

/// Parameters for the formal untweaked pSquare-hash permutation.
///
/// The construction consumes `t / 2` constants per Feistel round. This root-crate
/// implementation intentionally accepts only the two target instance shapes:
/// `t` in `{16, 24}` and `R = 52`.
#[derive(Clone, Debug)]
pub struct PSquareHashParams<F: FieldElement> {
    pub(crate) t: usize,
    pub(crate) rounds: usize,
    pub(crate) round_constants: Vec<F>,
}

impl<F: FieldElement> PSquareHashParams<F> {
    pub fn new(t: usize, rounds: usize, round_constants: &[F]) -> Self {
        assert!(
            SUPPORTED_WIDTHS.contains(&t),
            "pSquare-hash width must be 16 or 24"
        );
        assert_eq!(
            rounds, PSQUAREHASH_ROUNDS,
            "pSquare-hash formal instances use exactly 52 rounds"
        );
        assert_eq!(
            round_constants.len(),
            t / 2 * rounds,
            "pSquare-hash needs t / 2 round constants per round"
        );

        Self {
            t,
            rounds,
            round_constants: round_constants.to_vec(),
        }
    }
}

/// The formal untweaked pSquare-hash permutation.
#[derive(Clone, Debug)]
pub struct PSquareHash<F: FieldElement> {
    params: Arc<PSquareHashParams<F>>,
}

impl<F: FieldElement> PSquareHash<F> {
    pub fn new(params: &Arc<PSquareHashParams<F>>) -> Self {
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

    /// Allocating wrapper around the double-round implementation.
    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let mut state = input.to_vec();
        self.permutation_in_place(&mut state);
        state
    }

    /// Evaluate two consecutive Feistel rounds at a time.
    ///
    /// Pairing the rounds avoids materializing the half-state rotations performed
    /// by the literal reference algorithm.
    pub fn permutation_in_place(&self, state: &mut [F]) {
        self.validate_state(state);
        apply_matrix_in_place(state);

        for double_round in 0..self.params.rounds / 2 {
            let start = double_round * self.params.t;
            double_feistel_round(
                state,
                &self.params.round_constants[start..start + self.params.t],
            );
        }

        apply_matrix_in_place(state);
    }

    /// Literal one-Feistel-round-at-a-time implementation used as an audit oracle.
    pub fn permutation_reference(&self, input: &[F]) -> Vec<F> {
        self.validate_state(input);
        let constants_per_round = self.params.t / 2;
        let mut state = apply_matrix(input);

        for round in 0..self.params.rounds {
            let start = round * constants_per_round;
            state = feistel_round_reference(
                &state,
                &self.params.round_constants[start..start + constants_per_round],
            );
        }

        apply_matrix(&state)
    }

    fn validate_state(&self, state: &[F]) {
        assert_eq!(
            state.len(),
            self.params.t,
            "pSquare-hash input width does not match its parameters"
        );
    }
}

#[inline(always)]
fn apply_matrix<F: FieldElement>(state: &[F]) -> Vec<F> {
    let t = state.len();
    let half = t / 2;
    let mut result = vec![F::zero(); t];

    for i in 0..half {
        result[i] = state[i].clone();
        result[i].add_assign(&state[half + i]);

        result[half + i] = state[i].clone();
        result[half + i].double();
        result[half + i].add_assign(&state[half + i]);
    }

    result
}

#[inline(always)]
fn apply_matrix_in_place<F: FieldElement>(state: &mut [F]) {
    let half = state.len() / 2;

    for i in 0..half {
        let left = state[i].clone();
        let right = state[half + i].clone();

        state[i].add_assign(&right);

        let mut lower = left;
        lower.double();
        lower.add_assign(&right);
        state[half + i] = lower;
    }
}

#[inline(always)]
fn non_linear_layer<F: FieldElement>(state: &[F], constants: &[F]) -> [F; 2] {
    debug_assert_eq!(state.len(), 2);
    debug_assert_eq!(constants.len(), 2);

    let mut y_1 = state[1].clone();
    y_1.add_assign(&constants[0]);

    let mut y_3 = y_1.clone();
    y_3.add_assign(&constants[1]);

    y_1.square();
    let mut y_2 = state[0].clone();
    y_2.add_assign(&y_1);

    y_3.add_assign(&y_2);
    let mut y_5 = y_3.clone();

    y_3.square();
    let mut y_4 = y_2;
    y_4.add_assign(&y_3);
    y_5.add_assign(&y_4);

    [y_4, y_5]
}

#[inline(always)]
fn linear_combination<F: FieldElement>(state: &[F]) -> [F; 2] {
    debug_assert_eq!(state.len(), 2);

    let mut z_1 = state[0].clone();
    z_1.add_assign(&state[1]);
    let mut z_0 = z_1.clone();
    z_0.add_assign(&state[0]);

    [z_0, z_1]
}

fn feistel_round_reference<F: FieldElement>(state: &[F], constants: &[F]) -> Vec<F> {
    let t = state.len();
    let half = t / 2;
    let pair_count = half / 2;
    debug_assert_eq!(constants.len(), half);

    let z = linear_combination(&state[half - 2..half]);
    let mut result = vec![F::zero(); t];

    for pair in 0..pair_count {
        let source = half - 2 * (pair + 1);
        let destination = 2 * pair;
        let y = non_linear_layer(
            &state[source..source + 2],
            &constants[destination..destination + 2],
        );

        for lane in 0..2 {
            let mut value = state[half + destination + lane].clone();
            value.add_assign(&y[lane]);

            if pair > 0 {
                value.add_assign(&z[lane]);
            }
            if pair + 1 == pair_count {
                for index in (2..half - 2).step_by(2) {
                    value.add_assign(&state[index + lane]);
                }
            }

            result[destination + lane] = value;
        }
    }

    result[half..].clone_from_slice(&state[..half]);
    result
}

#[inline(always)]
fn double_feistel_round<F: FieldElement>(state: &mut [F], constants: &[F]) {
    let half = state.len() / 2;
    debug_assert_eq!(constants.len(), state.len());

    feistel_half_in_place(state, 0, half, &constants[..half]);
    feistel_half_in_place(state, half, 0, &constants[half..]);
}

/// Update the destination half with one Feistel round while retaining the source
/// half in place. A second call with the halves reversed completes a double round.
#[inline(always)]
fn feistel_half_in_place<F: FieldElement>(
    state: &mut [F],
    source_offset: usize,
    destination_offset: usize,
    constants: &[F],
) {
    let half = state.len() / 2;
    let pair_count = half / 2;
    debug_assert_eq!(constants.len(), half);

    let z = linear_combination(&state[source_offset + half - 2..source_offset + half]);

    for pair in 0..pair_count {
        let source = source_offset + half - 2 * (pair + 1);
        let destination = destination_offset + 2 * pair;
        let constant_offset = 2 * pair;
        let y = non_linear_layer(
            &state[source..source + 2],
            &constants[constant_offset..constant_offset + 2],
        );

        for lane in 0..2 {
            let mut increment = y[lane].clone();

            if pair > 0 {
                increment.add_assign(&z[lane]);
            }
            if pair + 1 == pair_count {
                for index in (2..half - 2).step_by(2) {
                    increment.add_assign(&state[source_offset + index + lane]);
                }
            }

            state[destination + lane].add_assign(&increment);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{PSquareHash, PSquareHashParams, PSQUAREHASH_ROUNDS};
    use crate::fields::babybear::BabyBear;
    use crate::fields::koalabear::KoalaBear;
    use crate::fields::mersenne31::Mersenne31;
    use crate::fields::FieldElement;
    use crate::psquarehash::instances::{
        PSQUAREHASH_BABYBEAR_16_PARAMS, PSQUAREHASH_BABYBEAR_24_PARAMS,
        PSQUAREHASH_KOALABEAR_16_PARAMS, PSQUAREHASH_KOALABEAR_24_PARAMS,
        PSQUAREHASH_MERSENNE31_16_PARAMS, PSQUAREHASH_MERSENNE31_24_PARAMS,
    };
    use std::sync::Arc;

    fn m31_vec(values: &[u32]) -> Vec<Mersenne31> {
        values.iter().copied().map(Mersenne31::from_u32).collect()
    }

    fn deterministic_input<F: FieldElement>(t: usize, seed: u32) -> Vec<F> {
        let mut value = seed;
        (0..t)
            .map(|_| {
                value = value.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                F::from_u64(value as u64)
            })
            .collect()
    }

    fn assert_reference_matches_optimized<F: FieldElement>(
        permutation: &PSquareHash<F>,
        seed: u32,
    ) {
        let input = deterministic_input(permutation.get_t(), seed);
        let expected = permutation.permutation_reference(&input);
        let mut actual = input;
        permutation.permutation_in_place(&mut actual);
        assert_eq!(actual, expected);
    }

    fn assert_known_answer<F: FieldElement>(params: &Arc<PSquareHashParams<F>>, expected: &[u32]) {
        let permutation = PSquareHash::new(params);
        let input: Vec<F> = (1..=permutation.get_t())
            .map(|value| F::from_u64(value as u64))
            .collect();
        let expected: Vec<F> = expected
            .iter()
            .map(|&value| F::from_u64(value as u64))
            .collect();

        assert_eq!(permutation.permutation(&input), expected);
        assert_eq!(permutation.permutation_reference(&input), expected);
    }

    #[test]
    fn optimized_matches_reference_for_target_instances() {
        for seed in [0, 1, 42, 0xdead_beef] {
            assert_reference_matches_optimized(
                &PSquareHash::new(&PSQUAREHASH_MERSENNE31_16_PARAMS),
                seed,
            );
            assert_reference_matches_optimized(
                &PSquareHash::new(&PSQUAREHASH_MERSENNE31_24_PARAMS),
                seed,
            );
            assert_reference_matches_optimized(
                &PSquareHash::new(&PSQUAREHASH_KOALABEAR_16_PARAMS),
                seed,
            );
            assert_reference_matches_optimized(
                &PSquareHash::new(&PSQUAREHASH_KOALABEAR_24_PARAMS),
                seed,
            );
            assert_reference_matches_optimized(
                &PSquareHash::new(&PSQUAREHASH_BABYBEAR_16_PARAMS),
                seed,
            );
            assert_reference_matches_optimized(
                &PSquareHash::new(&PSQUAREHASH_BABYBEAR_24_PARAMS),
                seed,
            );
        }
    }

    #[test]
    fn babybear_and_koalabear_known_answers_match_sister_repository() {
        assert_known_answer::<BabyBear>(
            &PSQUAREHASH_BABYBEAR_16_PARAMS,
            &[
                0x1b6b0497, 0x6462ea8a, 0x70409d39, 0x6c78d9e2, 0x3d9e9492, 0x03c2e5e7, 0x2ca27f5f,
                0x54e31173, 0x2dcc8298, 0x07dd9cb1, 0x54f34043, 0x1397ee75, 0x66a2621a, 0x6f801c1c,
                0x2a6d28a7, 0x3b8e756c,
            ],
        );
        assert_known_answer::<BabyBear>(
            &PSQUAREHASH_BABYBEAR_24_PARAMS,
            &[
                0x3f707db3, 0x13caf186, 0x1dc3c5e2, 0x0c38e088, 0x0adf6caa, 0x3ac23390, 0x5d549cfb,
                0x652de625, 0x0039d8a3, 0x54b8cfeb, 0x3118f068, 0x680ce302, 0x2ce4f7f9, 0x3dd95f74,
                0x0d9fafe4, 0x0c4565ed, 0x6eb8f7f4, 0x35808336, 0x67751724, 0x66644061, 0x73df8340,
                0x086fab19, 0x71a1310c, 0x65c1546c,
            ],
        );
        assert_known_answer::<KoalaBear>(
            &PSQUAREHASH_KOALABEAR_16_PARAMS,
            &[
                0x5c217aaa, 0x0d0c6b5c, 0x02fa1e09, 0x1740b34b, 0x518476f7, 0x02be94c7, 0x65abdaed,
                0x3bcc7ad3, 0x18f7296c, 0x7d6cb5a9, 0x044e0b17, 0x71671829, 0x1c90a2d3, 0x62908852,
                0x1fcfdb92, 0x1032b8b2,
            ],
        );
        assert_known_answer::<KoalaBear>(
            &PSQUAREHASH_KOALABEAR_24_PARAMS,
            &[
                0x18b17681, 0x4e8cc133, 0x687df48b, 0x1f61518b, 0x2d937612, 0x7e202012, 0x309617f9,
                0x3a1bc942, 0x280e7bc3, 0x0af7ede9, 0x2a6545bf, 0x4ffab561, 0x2d634b54, 0x21a173f9,
                0x3b46a286, 0x28950375, 0x69b0435e, 0x7e0cad36, 0x11d6dc54, 0x4f6f6531, 0x3c16c4c3,
                0x59c5c3d5, 0x4149c1a9, 0x2f64902e,
            ],
        );
    }

    #[test]
    fn mersenne31_t16_known_answer() {
        let permutation = PSquareHash::new(&PSQUAREHASH_MERSENNE31_16_PARAMS);
        let input = m31_vec(&[
            0x78066d6b, 0x68a24eb4, 0x2d12aacd, 0x42bb7df4, 0x3a85ecf4, 0x010084b5, 0x28e3f4fb,
            0x41514a49, 0x0e904f42, 0x0981bfd9, 0x3309b9ac, 0x19f408ff, 0x1f3202d0, 0x2ebbcc8c,
            0x261b659f, 0x22171a32,
        ]);
        let expected = m31_vec(&[
            0x0d31650f, 0x5b324b40, 0x02fb8ac7, 0x555c9139, 0x00a60cba, 0x1b61b003, 0x33e1c0ad,
            0x48d970a2, 0x0876e39d, 0x3a6f9513, 0x36afab87, 0x4d85ef87, 0x277b1cee, 0x70debee3,
            0x1337395b, 0x35ea5bad,
        ]);

        assert_eq!(permutation.permutation(&input), expected);
        assert_eq!(permutation.permutation_reference(&input), expected);
    }

    #[test]
    fn mersenne31_t24_known_answer() {
        let permutation = PSquareHash::new(&PSQUAREHASH_MERSENNE31_24_PARAMS);
        let input = m31_vec(&[
            0x78066d6b, 0x68a24eb4, 0x2d12aacd, 0x42bb7df4, 0x3a85ecf4, 0x010084b5, 0x28e3f4fb,
            0x41514a49, 0x0e904f42, 0x0981bfd9, 0x3309b9ac, 0x19f408ff, 0x1f3202d0, 0x2ebbcc8c,
            0x261b659f, 0x22171a32, 0x2b77fbfb, 0x57d3e692, 0x47dbb2c4, 0x5f803d52, 0x7791f988,
            0x6988c314, 0x283918dd, 0x32a8ab7b,
        ]);
        let expected = m31_vec(&[
            0x2e50ad42, 0x117b015c, 0x0c4610fb, 0x636c99be, 0x3b2635cc, 0x15323786, 0x36ba41ac,
            0x788cb4d8, 0x0f7b751a, 0x7608c969, 0x6ff8eeda, 0x6ed27e30, 0x55e7b993, 0x63506d14,
            0x42032061, 0x31bbe5f5, 0x7179d761, 0x5871965a, 0x16497d76, 0x2237878d, 0x637ca7d2,
            0x13752294, 0x0831c440, 0x18bf0647,
        ]);

        assert_eq!(permutation.permutation(&input), expected);
        assert_eq!(permutation.permutation_reference(&input), expected);
    }

    #[test]
    #[should_panic(expected = "pSquare-hash width must be 16 or 24")]
    fn constructor_rejects_unsupported_width() {
        let constants = vec![Mersenne31::zero(); 8 / 2 * PSQUAREHASH_ROUNDS];
        let _ = PSquareHashParams::new(8, PSQUAREHASH_ROUNDS, &constants);
    }

    #[test]
    #[should_panic(expected = "pSquare-hash formal instances use exactly 52 rounds")]
    fn constructor_rejects_wrong_round_count() {
        let constants = vec![Mersenne31::zero(); 16 / 2 * 51];
        let _ = PSquareHashParams::new(16, 51, &constants);
    }

    #[test]
    #[should_panic(expected = "pSquare-hash needs t / 2 round constants per round")]
    fn constructor_rejects_wrong_constant_count() {
        let constants = vec![Mersenne31::zero(); 16 / 2 * PSQUAREHASH_ROUNDS - 1];
        let _ = PSquareHashParams::new(16, PSQUAREHASH_ROUNDS, &constants);
    }
}
