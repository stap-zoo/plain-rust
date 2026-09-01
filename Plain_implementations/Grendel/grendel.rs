//! Plain Grendel permutation over prime fields.
//!
//! The paper specifies how to derive parameters but does not publish a recommended
//! instance table. The exported instances in this repository are reference-derived
//! at `kappa = 128`, including their round counts and the paper's 1.25 security margin.

use crate::fields::{FieldElement, PrimeField, PrimeFieldExt};
use crate::utils::modinv;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct GrendelParams<F: PrimeField> {
    pub(crate) t: usize,
    pub(crate) alpha: u64,
    pub(crate) rounds: usize,
    pub(crate) mds: Vec<Vec<F>>,
    pub(crate) mds_inv: Vec<Vec<F>>,
    pub(crate) round_constants: Vec<Vec<F>>,
    pub(crate) legendre_exp: Vec<u64>,
    pub(crate) sbox_inv_exp: Vec<u64>,
}

impl<F: PrimeField> GrendelParams<F> {
    pub fn new(
        t: usize,
        alpha: u64,
        rounds: usize,
        mds: &[Vec<F>],
        round_constants: &[Vec<F>],
    ) -> Self {
        assert!(t >= 2, "Grendel state width must be at least two");
        assert!(rounds >= 1, "Grendel must have at least one round");
        assert_eq!(mds.len(), t, "Grendel MDS matrix must have t rows");
        for row in mds {
            assert_eq!(row.len(), t, "Grendel MDS matrix must be t by t");
        }
        assert_eq!(
            round_constants.len(),
            rounds,
            "Grendel needs one constant vector per round"
        );
        for constants in round_constants {
            assert_eq!(
                constants.len(),
                t,
                "each Grendel constant vector must have width t"
            );
        }

        let modulus_minus_one = F::modulus() - BigUint::one();
        let legendre: BigUint = &modulus_minus_one >> 1usize;
        let sbox_exp = &legendre + BigUint::from(alpha);
        assert_eq!(
            gcd(sbox_exp.clone(), modulus_minus_one.clone()),
            BigUint::one(),
            "Grendel S-box exponent must be invertible modulo p - 1"
        );
        let sbox_inv = modinv(&sbox_exp, &modulus_minus_one);

        Self {
            t,
            alpha,
            rounds,
            mds: mds.to_vec(),
            mds_inv: invert_matrix(mds),
            round_constants: round_constants.to_vec(),
            legendre_exp: legendre.to_u64_digits(),
            sbox_inv_exp: sbox_inv.to_u64_digits(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Grendel<F: PrimeFieldExt> {
    params: Arc<GrendelParams<F>>,
}

impl<F: PrimeFieldExt> Grendel<F> {
    pub fn new(params: &Arc<GrendelParams<F>>) -> Self {
        Self {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn get_alpha(&self) -> u64 {
        self.params.alpha
    }

    pub fn get_rounds(&self) -> usize {
        self.params.rounds
    }

    /// Native permutation path, using the Jacobi algorithm for the Legendre symbol.
    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        self.validate_state(input);
        let mut state = input.to_vec();

        for round in 0..self.params.rounds {
            for value in &mut state {
                *value = self.sbox_fast(value);
            }
            state = matvec(&self.params.mds, &state);
            add_constants(&mut state, &self.params.round_constants[round]);
        }
        state
    }

    /// Literal reference path: compute `chi(x)` through Euler's criterion.
    pub fn permutation_reference(&self, input: &[F]) -> Vec<F> {
        self.validate_state(input);
        let mut state = input.to_vec();

        for round in 0..self.params.rounds {
            for value in &mut state {
                *value = self.sbox_reference(value);
            }
            state = matvec(&self.params.mds, &state);
            add_constants(&mut state, &self.params.round_constants[round]);
        }
        state
    }

    pub fn permutation_inv(&self, input: &[F]) -> Vec<F> {
        self.validate_state(input);
        let mut state = input.to_vec();

        for round in (0..self.params.rounds).rev() {
            sub_constants(&mut state, &self.params.round_constants[round]);
            state = matvec(&self.params.mds_inv, &state);
            for value in &mut state {
                *value = value.pow_words_le(&self.params.sbox_inv_exp);
            }
        }
        state
    }

    fn validate_state(&self, state: &[F]) {
        assert_eq!(
            state.len(),
            self.params.t,
            "Grendel input width does not match its parameters"
        );
    }

    fn sbox_fast(&self, value: &F) -> F {
        let powered = value.pow_u64(self.params.alpha);
        if jacobi_symbol(value.to_biguint(), F::modulus()) < 0 {
            powered.negate()
        } else {
            powered
        }
    }

    fn sbox_reference(&self, value: &F) -> F {
        let mut result = value.pow_u64(self.params.alpha);
        result.mul_assign(&value.pow_words_le(&self.params.legendre_exp));
        result
    }
}

fn matvec<F: FieldElement>(matrix: &[Vec<F>], vector: &[F]) -> Vec<F> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(vector)
                .fold(F::zero(), |mut acc, (coefficient, value)| {
                    let mut term = coefficient.clone();
                    term.mul_assign(value);
                    acc.add_assign(&term);
                    acc
                })
        })
        .collect()
}

fn add_constants<F: FieldElement>(state: &mut [F], constants: &[F]) {
    for (value, constant) in state.iter_mut().zip(constants) {
        value.add_assign(constant);
    }
}

fn sub_constants<F: FieldElement>(state: &mut [F], constants: &[F]) {
    for (value, constant) in state.iter_mut().zip(constants) {
        value.sub_assign(constant);
    }
}

fn invert_matrix<F: PrimeField>(matrix: &[Vec<F>]) -> Vec<Vec<F>> {
    let t = matrix.len();
    let mut augmented = vec![vec![F::zero(); 2 * t]; t];
    for i in 0..t {
        augmented[i][..t].clone_from_slice(&matrix[i]);
        augmented[i][t + i] = F::one();
    }

    let inverse_exp = (F::modulus() - BigUint::from(2u64)).to_u64_digits();
    for col in 0..t {
        let pivot = (col..t)
            .find(|&row| augmented[row][col] != F::zero())
            .expect("Grendel MDS matrix must be invertible");
        augmented.swap(col, pivot);

        let pivot_inv = augmented[col][col].pow_words_le(&inverse_exp);
        for entry in &mut augmented[col] {
            entry.mul_assign(&pivot_inv);
        }

        for row in 0..t {
            if row == col {
                continue;
            }
            let factor = augmented[row][col].clone();
            if factor == F::zero() {
                continue;
            }
            for j in 0..2 * t {
                let mut term = augmented[col][j].clone();
                term.mul_assign(&factor);
                augmented[row][j].sub_assign(&term);
            }
        }
    }

    augmented.into_iter().map(|row| row[t..].to_vec()).collect()
}

fn gcd(mut left: BigUint, mut right: BigUint) -> BigUint {
    while !right.is_zero() {
        let remainder = left % &right;
        left = right;
        right = remainder;
    }
    left
}

/// Binary Jacobi symbol `(a / n)` for positive odd `n`.
fn jacobi_symbol(mut a: BigUint, mut n: BigUint) -> i8 {
    assert!(
        !n.is_zero() && n.bit(0),
        "Jacobi modulus must be positive and odd"
    );
    a %= &n;
    let mut sign = 1i8;

    while !a.is_zero() {
        while !a.bit(0) {
            a >>= 1;
            let n_mod_8 = (&n % BigUint::from(8u64)).to_u64_digits()[0];
            if n_mod_8 == 3 || n_mod_8 == 5 {
                sign = -sign;
            }
        }

        std::mem::swap(&mut a, &mut n);
        if (&a % BigUint::from(4u64)) == BigUint::from(3u64)
            && (&n % BigUint::from(4u64)) == BigUint::from(3u64)
        {
            sign = -sign;
        }
        a %= &n;
    }

    if n == BigUint::one() {
        sign
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::bls12_381::Bls12_381;
    use crate::fields::goldilocks::Goldilocks;
    use crate::grendel::instances::*;

    fn field_vec<F: PrimeField>(values: &[&str]) -> Vec<F> {
        values
            .iter()
            .map(|value| {
                F::from_biguint(
                    &BigUint::parse_bytes(value.as_bytes(), 10)
                        .expect("KAT field element must be decimal"),
                )
            })
            .collect()
    }

    fn check_kat<F: PrimeFieldExt>(params: &Arc<GrendelParams<F>>, expected: &[&str]) {
        let permutation = Grendel::new(params);
        let input = (1..=permutation.get_t())
            .map(|value| F::from_u64(value as u64))
            .collect::<Vec<_>>();
        let expected = field_vec::<F>(expected);

        assert_eq!(permutation.permutation(&input), expected);
        assert_eq!(permutation.permutation_reference(&input), expected);
    }

    #[test]
    fn reference_known_answer_vectors() {
        check_kat(
            &GRENDEL_BN254_3_PARAMS,
            &[
                "15995791472288294308174169687998331711129185833179527498191054822616746540318",
                "1277715001056385596364846198202477548633516393666932412963108054529318007059",
                "18340502957806533494344771706393090410841925798215504224892478240497454278641",
            ],
        );
        check_kat(
            &GRENDEL_BLS12_381_3_PARAMS,
            &[
                "9434234317712565218121184097139083820371497410578881151920995205442822400521",
                "24298086792365008571996165896764301399721420572426783449784437920605866399811",
                "49918702889014010387344552875446168901658341302599961293432898798414336052092",
            ],
        );
        check_kat(
            &GRENDEL_BN254_4_CURRENT_PARAMS,
            &[
                "13206260123643548167057821422681843389124331583678963099474887108000716342338",
                "15945999528250322765225831664723251869126099531593449439914855171711522098528",
                "15637910753041344178139080674759404851635589135818682821907752110249864303732",
                "17872999515418434829882019995023785167607565254232172343013161080602228187619",
            ],
        );
        check_kat(
            &GRENDEL_BLS12_381_4_CURRENT_PARAMS,
            &[
                "20820025499377825222754586235024540014569409009297806375598851919555741006230",
                "20554745102470582060180954574648489231967250807396942781777266093644510996156",
                "30630179309470311194352383575991366828813121001168205609925548715464340371081",
                "37457310864457231113051423226476773373012903179701527760519589254119386953168",
            ],
        );
        check_kat(
            &GRENDEL_GOLDILOCKS_8_PARAMS,
            &[
                "1773155256352791889",
                "9608147794966071694",
                "6936904078599778261",
                "7924792274788895277",
                "5260515945985951191",
                "7332056128811092845",
                "9085930165439481026",
                "6602104553378359669",
            ],
        );
        check_kat(
            &GRENDEL_GOLDILOCKS_12_PARAMS,
            &[
                "5035495851539612286",
                "14631084273934423798",
                "18146424514498527994",
                "505625766681736420",
                "8577393251852828811",
                "14975937896620890716",
                "10378365540883483745",
                "3356291182947947242",
                "3469645843760651611",
                "16053439067072341888",
                "794727091748218690",
                "3837276950359648869",
            ],
        );
        check_kat(
            &GRENDEL_MERSENNE31_16_PARAMS,
            &[
                "1548498469",
                "220518536",
                "758626964",
                "577542462",
                "1574627403",
                "276282669",
                "1862908292",
                "1933812200",
                "1625682337",
                "113479202",
                "1940706437",
                "165255642",
                "1248558055",
                "1689177701",
                "1971099579",
                "488955104",
            ],
        );
        check_kat(
            &GRENDEL_MERSENNE31_24_PARAMS,
            &[
                "1280816985",
                "1351205686",
                "1528332845",
                "326378604",
                "1471674389",
                "1293292128",
                "1887745967",
                "287816800",
                "873489363",
                "900633908",
                "1794157424",
                "2036545510",
                "878844316",
                "1324930115",
                "1122806003",
                "944629631",
                "1759198900",
                "406598449",
                "475434585",
                "435411831",
                "1768564338",
                "1880126791",
                "1963007527",
                "1018146772",
            ],
        );
        check_kat(
            &GRENDEL_BABYBEAR_16_PARAMS,
            &[
                "240603057",
                "141589259",
                "1880134246",
                "1665694529",
                "1124922814",
                "1669103594",
                "884382146",
                "1341318542",
                "1342673176",
                "1160872775",
                "1397752431",
                "1753398632",
                "1920086178",
                "818494683",
                "1235172711",
                "635114624",
            ],
        );
        check_kat(
            &GRENDEL_BABYBEAR_24_PARAMS,
            &[
                "137275462",
                "1374467081",
                "1430162970",
                "230430939",
                "1104549180",
                "1506796573",
                "937619567",
                "460451141",
                "205405579",
                "1225067243",
                "1693583162",
                "788967721",
                "872445550",
                "292788677",
                "916623559",
                "1198437579",
                "546755411",
                "1163784844",
                "1167043629",
                "1264685387",
                "947258984",
                "1593783011",
                "1364311452",
                "700394408",
            ],
        );
        check_kat(
            &GRENDEL_KOALABEAR_16_PARAMS,
            &[
                "1263694629",
                "1576844970",
                "1895926089",
                "1939208923",
                "14701756",
                "970061525",
                "967724396",
                "882246125",
                "345122034",
                "1697437113",
                "1339074409",
                "1581380265",
                "1430790338",
                "987239116",
                "632743990",
                "734764946",
            ],
        );
        check_kat(
            &GRENDEL_KOALABEAR_24_PARAMS,
            &[
                "1605452543",
                "1828409724",
                "1260515453",
                "413516661",
                "1719716141",
                "413530239",
                "9359488",
                "104201057",
                "786670888",
                "1816206251",
                "235326317",
                "1450430187",
                "1900318501",
                "601343489",
                "1678901715",
                "1489078892",
                "1741901188",
                "400046932",
                "1113778059",
                "828299756",
                "1836365924",
                "748367581",
                "1254195791",
                "1595321627",
            ],
        );
    }

    fn check_paths_and_inverse<F: PrimeFieldExt>(params: &Arc<GrendelParams<F>>) {
        let permutation = Grendel::new(params);
        let mut lcg = 0x9e37_79b9_7f4a_7c15u64;
        let input = (0..permutation.get_t())
            .map(|_| {
                lcg = lcg
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                F::from_u64(lcg)
            })
            .collect::<Vec<_>>();

        let fast = permutation.permutation(&input);
        assert_eq!(fast, permutation.permutation_reference(&input));
        assert_eq!(permutation.permutation_inv(&fast), input);
    }

    #[test]
    fn jacobi_and_reference_paths_agree_and_round_trip() {
        check_paths_and_inverse(&GRENDEL_BN254_3_PARAMS);
        check_paths_and_inverse(&GRENDEL_BLS12_381_3_PARAMS);
        check_paths_and_inverse(&GRENDEL_BN254_4_CURRENT_PARAMS);
        check_paths_and_inverse(&GRENDEL_BN254_4_ORIGINAL_PARAMS);
        check_paths_and_inverse(&GRENDEL_BLS12_381_4_CURRENT_PARAMS);
        check_paths_and_inverse(&GRENDEL_BLS12_381_4_ORIGINAL_PARAMS);
        check_paths_and_inverse(&GRENDEL_GOLDILOCKS_8_PARAMS);
        check_paths_and_inverse(&GRENDEL_GOLDILOCKS_12_PARAMS);
        check_paths_and_inverse(&GRENDEL_MERSENNE31_16_PARAMS);
        check_paths_and_inverse(&GRENDEL_MERSENNE31_24_PARAMS);
        check_paths_and_inverse(&GRENDEL_BABYBEAR_16_PARAMS);
        check_paths_and_inverse(&GRENDEL_BABYBEAR_24_PARAMS);
        check_paths_and_inverse(&GRENDEL_KOALABEAR_16_PARAMS);
        check_paths_and_inverse(&GRENDEL_KOALABEAR_24_PARAMS);
    }

    #[test]
    fn zero_is_fixed_by_the_sbox() {
        let permutation = Grendel::new(&GRENDEL_GOLDILOCKS_8_PARAMS);
        let zero = Goldilocks::zero();
        assert_eq!(permutation.sbox_fast(&zero), zero);
        assert_eq!(permutation.sbox_reference(&zero), zero);
    }

    #[test]
    fn derived_alpha_and_round_table_is_pinned() {
        assert_eq!(
            (GRENDEL_BN254_3_PARAMS.alpha, GRENDEL_BN254_3_PARAMS.rounds),
            (5, 24)
        );
        assert_eq!(
            (
                GRENDEL_BLS12_381_3_PARAMS.alpha,
                GRENDEL_BLS12_381_3_PARAMS.rounds
            ),
            (5, 24)
        );
        assert_eq!(
            (
                GRENDEL_BN254_4_CURRENT_PARAMS.alpha,
                GRENDEL_BN254_4_CURRENT_PARAMS.rounds
            ),
            (5, 23)
        );
        assert_eq!(
            (
                GRENDEL_BN254_4_ORIGINAL_PARAMS.alpha,
                GRENDEL_BN254_4_ORIGINAL_PARAMS.rounds
            ),
            (5, 12)
        );
        assert_eq!(
            (
                GRENDEL_BLS12_381_4_CURRENT_PARAMS.alpha,
                GRENDEL_BLS12_381_4_CURRENT_PARAMS.rounds
            ),
            (5, 23)
        );
        assert_eq!(
            (
                GRENDEL_BLS12_381_4_ORIGINAL_PARAMS.alpha,
                GRENDEL_BLS12_381_4_ORIGINAL_PARAMS.rounds
            ),
            (5, 12)
        );
        assert_eq!(
            (
                GRENDEL_GOLDILOCKS_8_PARAMS.alpha,
                GRENDEL_GOLDILOCKS_8_PARAMS.rounds
            ),
            (7, 15)
        );
        assert_eq!(
            (
                GRENDEL_GOLDILOCKS_12_PARAMS.alpha,
                GRENDEL_GOLDILOCKS_12_PARAMS.rounds
            ),
            (7, 11)
        );
        assert_eq!(
            (
                GRENDEL_MERSENNE31_16_PARAMS.alpha,
                GRENDEL_MERSENNE31_16_PARAMS.rounds
            ),
            (2, 12)
        );
        assert_eq!(
            (
                GRENDEL_MERSENNE31_24_PARAMS.alpha,
                GRENDEL_MERSENNE31_24_PARAMS.rounds
            ),
            (2, 12)
        );
        assert_eq!(
            (
                GRENDEL_BABYBEAR_16_PARAMS.alpha,
                GRENDEL_BABYBEAR_16_PARAMS.rounds
            ),
            (7, 14)
        );
        assert_eq!(
            (
                GRENDEL_BABYBEAR_24_PARAMS.alpha,
                GRENDEL_BABYBEAR_24_PARAMS.rounds
            ),
            (7, 14)
        );
        assert_eq!(
            (
                GRENDEL_KOALABEAR_16_PARAMS.alpha,
                GRENDEL_KOALABEAR_16_PARAMS.rounds
            ),
            (3, 12)
        );
        assert_eq!(
            (
                GRENDEL_KOALABEAR_24_PARAMS.alpha,
                GRENDEL_KOALABEAR_24_PARAMS.rounds
            ),
            (3, 12)
        );
    }

    fn assert_same_non_round_profile<F: PrimeField>(
        current: &GrendelParams<F>,
        original: &GrendelParams<F>,
    ) {
        assert_eq!(current.t, original.t);
        assert_eq!(current.alpha, original.alpha);
        assert_eq!(current.mds, original.mds);
        assert_eq!(current.mds_inv, original.mds_inv);
        assert_eq!(
            &current.round_constants[..original.rounds],
            original.round_constants
        );
    }

    #[test]
    fn paired_profiles_share_all_non_round_choices() {
        assert_same_non_round_profile(
            &GRENDEL_BN254_4_CURRENT_PARAMS,
            &GRENDEL_BN254_4_ORIGINAL_PARAMS,
        );
        assert_same_non_round_profile(
            &GRENDEL_BLS12_381_4_CURRENT_PARAMS,
            &GRENDEL_BLS12_381_4_ORIGINAL_PARAMS,
        );
    }

    #[test]
    fn bls12_381_mds_matches_the_reference_export() {
        let expected = field_vec::<Bls12_381>(&[
            "52435875175126190479447740508185965837690552500527637822603658699938581066864",
            "137200",
            "52435875175126190479447740508185965837690552500527637822603658699938581164563",
            "400",
            "52435875175126190479447740508185965837690552500527637822603658699938534124913",
            "54762351",
            "52435875175126190479447740508185965837690552500527637822603658699938573341713",
            "140050",
            "52435875175126190479447740508185965837690552500527637822603658699922104442063",
            "19167800400",
            "52435875175126190479447740508185965837690552500527637822603658699935841949364",
            "48177200",
            "52435875175126190479447740508185965837690552500527637822603658694270581781713",
            "6593435097550",
            "52435875175126190479447740508185965837690552500527637822603658698996613844913",
            "16531644851",
        ]);
        let actual = GRENDEL_BLS12_381_4_CURRENT_PARAMS
            .mds
            .iter()
            .flatten()
            .cloned()
            .collect::<Vec<_>>();
        assert_eq!(actual, expected);
    }

    fn valid_t2_params() -> (Vec<Vec<Goldilocks>>, Vec<Vec<Goldilocks>>) {
        (
            vec![
                vec![Goldilocks::one(), Goldilocks::zero()],
                vec![Goldilocks::zero(), Goldilocks::one()],
            ],
            vec![vec![Goldilocks::zero(), Goldilocks::zero()]],
        )
    }

    #[test]
    #[should_panic(expected = "at least two")]
    fn rejects_width_below_two() {
        GrendelParams::<Goldilocks>::new(
            1,
            7,
            1,
            &[vec![Goldilocks::one()]],
            &[vec![Goldilocks::zero()]],
        );
    }

    #[test]
    #[should_panic(expected = "t by t")]
    fn rejects_wrong_matrix_shape() {
        let (_, constants) = valid_t2_params();
        GrendelParams::<Goldilocks>::new(
            2,
            7,
            1,
            &[vec![Goldilocks::one()], vec![Goldilocks::zero()]],
            &constants,
        );
    }

    #[test]
    #[should_panic(expected = "constant vector")]
    fn rejects_wrong_constant_shape() {
        let (matrix, _) = valid_t2_params();
        GrendelParams::<Goldilocks>::new(2, 7, 1, &matrix, &[vec![Goldilocks::zero()]]);
    }
}
