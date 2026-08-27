use super::polocolo::PolocoloParams;
use super::polocolo_luts::{
    POLOCOLO_BLS12_381_T3_M1024_LUT_WORDS, POLOCOLO_BN254_T3_M1024_LUT_WORDS,
};
use super::polocolo_sigma512::POLOCOLO_SIGMA_512;
use crate::fields::bls12_381::Bls12_381;
use crate::fields::bn254::Bn254;
use crate::fields::{
    biguint_from_limbs_le, biguint_to_limbs_le_4, FieldElement, PrimeField, PrimeFieldWords,
};
use lazy_static::lazy_static;
use num_bigint::BigUint;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake128;
use std::collections::HashMap;
use std::sync::Arc;

const T: usize = 3;
const ROUNDS: usize = 6;
const M: usize = 1024;

fn f_from_words<F: PrimeField>(words: [u64; 4]) -> F {
    F::from_biguint(&biguint_from_limbs_le(&words))
}

lazy_static! {
    pub static ref POLOCOLO_BN254_4_PARAMS: Arc<PolocoloParams<Bn254>> =
        derive_t4_params(5, "BN254");
    pub static ref POLOCOLO_BLS12_381_4_PARAMS: Arc<PolocoloParams<Bls12_381>> =
        derive_t4_params(7, "BLS12");
}

fn mds_t3<F: FieldElement>() -> Vec<Vec<F>> {
    vec![
        vec![F::from_u64(2), F::from_u64(1), F::from_u64(1)],
        vec![F::from_u64(1), F::from_u64(2), F::from_u64(1)],
        vec![F::from_u64(1), F::from_u64(1), F::from_u64(2)],
    ]
}

fn mds_t4<F: FieldElement>() -> Vec<Vec<F>> {
    [[5u64, 7, 1, 3], [4, 6, 1, 1], [1, 3, 5, 7], [1, 1, 4, 6]]
        .iter()
        .map(|row| row.iter().map(|value| F::from_u64(*value)).collect())
        .collect()
}

fn sample_bitshift_be<F: PrimeField>(reader: &mut dyn XofReader) -> F {
    let modulus = F::modulus();
    let bits = modulus.bits() as usize;
    let width = bits.div_ceil(8);
    let shift = (8 - bits % 8) % 8;
    loop {
        let mut bytes = vec![0u8; width];
        reader.read(&mut bytes);
        bytes[0] >>= shift;
        let candidate = BigUint::from_bytes_be(&bytes);
        if candidate < modulus {
            return F::from_biguint(&candidate);
        }
    }
}

fn derive_t4_params<F: PrimeField + PrimeFieldWords>(
    generator: u64,
    field_label: &str,
) -> Arc<PolocoloParams<F>> {
    const T4: usize = 4;
    const M4: usize = 512;
    const ROUNDS4: usize = 5;

    let modulus = F::modulus();
    let ann_big = (modulus.clone() - BigUint::from(1u8)) / BigUint::from(M4);
    let ann = biguint_to_limbs_le_4(&ann_big);
    let modulus_minus_two = biguint_to_limbs_le_4(&(modulus - BigUint::from(2u8)));
    let mds = mds_t4::<F>();

    let seed = format!("Polocolo-{M4}-{ROUNDS4}-{T4}-{field_label}");
    let mut shake = Shake128::default();
    shake.update(seed.as_bytes());
    let mut reader = shake.finalize_xof();
    let round_constants = (0..ROUNDS4)
        .map(|_| {
            (0..T4)
                .map(|_| sample_bitshift_be::<F>(&mut reader))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let g = F::from_u64(generator);
    let gk = g.pow_words_le(&ann);
    let mut residue = F::one();
    let mut lut = HashMap::with_capacity(M4 + 1);
    lut.insert(F::zero().to_words_le(), F::zero());
    for (r, sigma) in POLOCOLO_SIGMA_512.iter().enumerate() {
        let exponent = ((M4 + 1) * r + sigma) as u64;
        lut.insert(residue.to_words_le(), g.pow_u64(exponent));
        residue.mul_assign(&gk);
    }

    Arc::new(PolocoloParams::new(
        T4,
        M4,
        ann,
        ROUNDS4,
        &mds,
        &round_constants,
        &lut,
        modulus_minus_two,
    ))
}

fn build_lut_and_modulus_minus_two<F: PrimeField>(
    lut_words: &[([u64; 4], [u64; 4])],
) -> (HashMap<[u64; 4], F>, [u64; 4]) {
    let mut lut = HashMap::with_capacity(lut_words.len());
    for (lut_in_words, lut_out_words) in lut_words.iter() {
        lut.insert(*lut_in_words, f_from_words::<F>(*lut_out_words));
    }

    let modulus = F::modulus();
    let modulus_minus_two = &modulus - BigUint::from(2u32);
    (lut, biguint_to_limbs_le_4(&modulus_minus_two))
}

lazy_static! {
    pub static ref POLOCOLO_BN254_3_PARAMS: Arc<PolocoloParams<Bn254>> = {
        type Scalar = Bn254;

        let ann = [
            0x2450f87d64fc0000u64,
            0x16ca0cfa121e6e5cu64,
            0x0a6e14116da06056u64,
            0x000c19139cb84c68u64,
        ];
        let mds = mds_t3::<Scalar>();
        let round_constants: Vec<Vec<Scalar>> = vec![
            vec![
                f_from_words::<Scalar>([
                    0xdfaaa8988b5d4763u64,
                    0x78eff157836b9cb7u64,
                    0x9d7998b55532c1c7u64,
                    0x21dd996b4a62c949u64,
                ]),
                f_from_words::<Scalar>([
                    0xf066fc429beac1c5u64,
                    0xd28b157416209c0eu64,
                    0x9b4d67c952ed87a7u64,
                    0x144228e01e23a698u64,
                ]),
                f_from_words::<Scalar>([
                    0x8ce8a0751fb1b2d5u64,
                    0xa8bd7502ec47fe56u64,
                    0xf640b91f453aa47fu64,
                    0x0631cbadc9fadd06u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0xd712e87fca7e78b7u64,
                    0x40a8d3744998ac4au64,
                    0x44665070b48a570du64,
                    0x2936aaa67ca846a5u64,
                ]),
                f_from_words::<Scalar>([
                    0x52c7424c38655d74u64,
                    0x438f60dececd5fecu64,
                    0x1af5a93354deb08bu64,
                    0x14653f268622bdb9u64,
                ]),
                f_from_words::<Scalar>([
                    0x2616209d8d696093u64,
                    0xe8ed70eaeb43a109u64,
                    0x3a2b70983c50c968u64,
                    0x0af70780f6be79f1u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0xe4cc64579326dc62u64,
                    0x03cec96f4cac1696u64,
                    0x4bcc5cf75b339c75u64,
                    0x2b7e535e5c81fa75u64,
                ]),
                f_from_words::<Scalar>([
                    0xa86e3c2a21f3f0b7u64,
                    0xd09882829e2a3487u64,
                    0x88ee6fdc19210c9cu64,
                    0x29a3b253beadcad3u64,
                ]),
                f_from_words::<Scalar>([
                    0xa25a80d8b60e017bu64,
                    0x2621a7ebf662207au64,
                    0xed5e3fe0652cde84u64,
                    0x23a38cee9f81e461u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0xbea045c9686dabbeu64,
                    0x41b6d4239c69b20cu64,
                    0xa4197a4436d3b1dcu64,
                    0x11fc5357634621ccu64,
                ]),
                f_from_words::<Scalar>([
                    0x3a3dad8fb66a96cfu64,
                    0xc47cc27bad98d239u64,
                    0x4402544463252665u64,
                    0x0fcf6e1e688ba7bau64,
                ]),
                f_from_words::<Scalar>([
                    0x2065203eeb6fa7e3u64,
                    0x41405a17b734e82cu64,
                    0x09fbc478ea04f747u64,
                    0x258f2e7f423b4648u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0x692768adce1a8709u64,
                    0x0ab804aaeac5ae77u64,
                    0xc47386c2457b6a39u64,
                    0x1a6dffb5af6ce61eu64,
                ]),
                f_from_words::<Scalar>([
                    0x630e0fb2e3e1cd5bu64,
                    0x9c44b608d4bbe503u64,
                    0x7a637dc1a0962864u64,
                    0x155658e4092c7658u64,
                ]),
                f_from_words::<Scalar>([
                    0x0b1c22343a4dbcf9u64,
                    0x338d7c1d14a1a817u64,
                    0xb4256b7bcfa9cb7au64,
                    0x10e89e2d08e17b68u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0xaac75bdda573a050u64,
                    0x4039f043b05f5f20u64,
                    0xd884b85a2692d65au64,
                    0x2ca2b58a89bd6e96u64,
                ]),
                f_from_words::<Scalar>([
                    0x5e79fbd6d5449b7cu64,
                    0x3bdde92a99d8e8ceu64,
                    0x487edfdd5aca8663u64,
                    0x25196201b4f47d2eu64,
                ]),
                f_from_words::<Scalar>([
                    0xb2ffd1868d01a080u64,
                    0xc1cdb0cce184eb7cu64,
                    0x47a4933588191a37u64,
                    0x0a4ac53202c4cfffu64,
                ]),
            ],
        ];

        let (lut, modulus_minus_two) =
            build_lut_and_modulus_minus_two::<Scalar>(&POLOCOLO_BN254_T3_M1024_LUT_WORDS);
        Arc::new(PolocoloParams::new(
            T,
            M,
            ann,
            ROUNDS,
            &mds,
            &round_constants,
            &lut,
            modulus_minus_two,
        ))
    };
    pub static ref POLOCOLO_BLS12_381_3_PARAMS: Arc<PolocoloParams<Bls12_381>> = {
        type Scalar = Bls12_381;

        let ann = [
            0xffbfffffffc00000u64,
            0x0154ef6900bfff96u64,
            0x520cce7602026876u64,
            0x001cfb69d4ca675fu64,
        ];
        let mds = mds_t3::<Scalar>();
        let round_constants: Vec<Vec<Scalar>> = vec![
            vec![
                f_from_words::<Scalar>([
                    0x80e8cad22b10f1a2u64,
                    0x03d3108651aa3928u64,
                    0xeebf0b080e950919u64,
                    0x09f0d43341118eecu64,
                ]),
                f_from_words::<Scalar>([
                    0x0f2c807de0941317u64,
                    0x57ebccccdb3c1df7u64,
                    0x8c8b0550bbee0382u64,
                    0x06220d2c73c88f8bu64,
                ]),
                f_from_words::<Scalar>([
                    0x489cd29f37caa730u64,
                    0x7a0c9de0f2818425u64,
                    0x7d99c737d28b127au64,
                    0x01a66aa665615458u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0x89133ec6d4eccdedu64,
                    0x516d3bc8a9b59391u64,
                    0x391bb7f6375123beu64,
                    0x166fe9de51094e7bu64,
                ]),
                f_from_words::<Scalar>([
                    0xa500cd6f6f57372du64,
                    0xf1b9909627580c62u64,
                    0xbda2ac3b9c388b6fu64,
                    0x598ecde013468f90u64,
                ]),
                f_from_words::<Scalar>([
                    0x39b88223192bfb01u64,
                    0x8e1f39f7b4a1a4d3u64,
                    0xb3f3d4e3d78bd07fu64,
                    0x63ad86431f7f297cu64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0x6dd5dc5f014cd033u64,
                    0x37ff4de9447e450du64,
                    0x3a7db2895f999bbeu64,
                    0x08010d2d51cc44a1u64,
                ]),
                f_from_words::<Scalar>([
                    0xef6f0359848c2805u64,
                    0x9160f4351d0a11cau64,
                    0xbb811bd5afcda80fu64,
                    0x08df5115cd3a91bau64,
                ]),
                f_from_words::<Scalar>([
                    0x6e91de04d244d48au64,
                    0xa36a493b910bb604u64,
                    0x3cb298b7cf0ae0b1u64,
                    0x2450546a0b2833f3u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0x335aebb2eb5ad6d9u64,
                    0x8e2245f68aea3c0eu64,
                    0xec739a5a7d9e6ceau64,
                    0x6aa680f4118007aeu64,
                ]),
                f_from_words::<Scalar>([
                    0xf8f0ff172303a90fu64,
                    0x0e93856514df5963u64,
                    0xd1faa999a5872bb6u64,
                    0x24e59d759fddb8d7u64,
                ]),
                f_from_words::<Scalar>([
                    0x373a60c4d59a05a0u64,
                    0x1f9ea637fb8fc8ceu64,
                    0xd39cd7f238fa2eaeu64,
                    0x2dccc660e814b318u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0x53c2663364e456edu64,
                    0xa6d0b69d5931fecfu64,
                    0x834159b570775222u64,
                    0x5153312c42e5dab3u64,
                ]),
                f_from_words::<Scalar>([
                    0x5a0f32e194c80c02u64,
                    0xd9d7fee43faf3914u64,
                    0x8495012a9f97766au64,
                    0x5a77e95babcd5215u64,
                ]),
                f_from_words::<Scalar>([
                    0xf8bd24dc766a58f5u64,
                    0xedac53e26716bdd9u64,
                    0x16726e41d139d24du64,
                    0x1a04f869d7349957u64,
                ]),
            ],
            vec![
                f_from_words::<Scalar>([
                    0xc4c3b9dd8175365cu64,
                    0xe84754f5f19cc162u64,
                    0x8c4357a129c724f9u64,
                    0x6078bd168ce702acu64,
                ]),
                f_from_words::<Scalar>([
                    0x08494b64da2b0e33u64,
                    0xbed90b51e7188e99u64,
                    0xb5d861469da7d91bu64,
                    0x2c9b0ff2c96360d4u64,
                ]),
                f_from_words::<Scalar>([
                    0x21c72345c592d810u64,
                    0xcbf6ae0277ef84a3u64,
                    0x798eeb74e2cced42u64,
                    0x2eff11a5d8aed46eu64,
                ]),
            ],
        ];

        let (lut, modulus_minus_two) =
            build_lut_and_modulus_minus_two::<Scalar>(&POLOCOLO_BLS12_381_T3_M1024_LUT_WORDS);
        Arc::new(PolocoloParams::new(
            T,
            M,
            ann,
            ROUNDS,
            &mds,
            &round_constants,
            &lut,
            modulus_minus_two,
        ))
    };
}
