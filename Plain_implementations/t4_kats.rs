use crate::anemoi::anemoi::Anemoi;
use crate::anemoi::instances::ANEMOI_BLS12_381_4_CURRENT_PARAMS;
use crate::arion::arion::Arion;
use crate::arion::instances::ARION_BLS12_381_4_PARAMS;
use crate::fields::bls12_381::Bls12_381;
use crate::fields::bn254::Bn254;
use crate::fields::FieldElement;
use crate::gmimc::gmimc::Gmimc;
use crate::gmimc::instances::GMIMC_BLS12_381_4_PARAMS;
use crate::gmimc2::gmimc2::Gmimc2;
use crate::gmimc2::instances::GMIMC2_BLS12_381_4_PARAMS;
use crate::grendel::grendel::Grendel;
use crate::grendel::instances::GRENDEL_BLS12_381_4_CURRENT_PARAMS;
use crate::griffin::griffin::Griffin;
use crate::griffin::instances::GRIFFIN_BLS12_381_4_CURRENT_PARAMS;
use crate::neptune::instances::NEPTUNE_BLS12_381_4_PARAMS;
use crate::neptune::neptune::Neptune;
use crate::polocolo::instances::POLOCOLO_BLS12_381_4_PARAMS;
use crate::polocolo::polocolo::Polocolo;
use crate::poseidon::instances::{POSEIDON_BLS12_381_4_PARAMS, POSEIDON_BN254_4_PARAMS};
use crate::poseidon::poseidon::Poseidon;
use crate::poseidon2::instances::{
    POSEIDON2_BLS12_381_4_PARAMS, POSEIDON2_BN254_4_GNARK_CRYPTO_PARAMS,
};
use crate::poseidon2::poseidon2::Poseidon2;
use crate::rescueprime::instances::RESCUE_PRIME_BLS12_381_4_CURRENT_PARAMS;
use crate::rescueprime::rescue_prime::RescuePrime;
use crate::skyscraper::instances::SKYSCRAPER_BLS12_381_2_PARAMS;
use crate::skyscraper::skyscraper::Skyscraper;
fn input() -> Vec<Bls12_381> {
    (1..=4).map(Bls12_381::from_u64).collect()
}
fn expected(values: &[&str]) -> Vec<Bls12_381> {
    values
        .iter()
        .map(|v| Bls12_381::from_hex(v).unwrap())
        .collect()
}

fn expected_bn254(values: &[&str]) -> Vec<Bn254> {
    values.iter().map(|v| Bn254::from_hex(v).unwrap()).collect()
}
#[test]
fn gmimc_bls12_t4_kat() {
    let permutation = Gmimc::new(&GMIMC_BLS12_381_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x5242950c941b9ef398645f2ca1ff0eec2a3353369640ef65d40e9cad1e1041e3",
            "0x540da6527b8d57ce9299210f287df0cf340b926515151d95622bddb79a67c4e0",
            "0x4b6345b31f7f7c5bb276a5e16efb55843e079703bd011fd43909f23317783a3a",
            "0x65c533e7ed793d61cd5f415d02040ea78b98561bf90d07c8b57591ac5e08f595"
        ])
    );
}
#[test]
fn gmimc2_bls12_t4_kat() {
    let permutation = Gmimc2::new(&GMIMC2_BLS12_381_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x5d0f952487bef007e440212648fe30f2a5fcf459d72d0d342084463f12d20d1e",
            "0x31e00229ba7605c5b00e7b16cc973cd7399a8c69065ad296939047d3e4efd0ec",
            "0x42c3a055ba3d42c66469f6e8ff0df54802b5218ca4d918701a8161df02625021",
            "0x461ed2e805b2272135ffbb93f2d9f6153aba1e7f7df5fe1f46ecdf20cc158c71"
        ])
    );
}
#[test]
fn poseidon_bls12_t4_kat() {
    let permutation = Poseidon::new(&POSEIDON_BLS12_381_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x6323c75327f7c80df5ae11235aa1185d6058bdab948a83240111a8ca6923dde4",
            "0x5b391ff979702978339267f5851f7fbdb4d712ffc12a54849d06252a2244426",
            "0x14c31cbb70cec2006b124b0fec4d9c0637281f958420203cc59f67dc7bad53f7",
            "0x1a4e3dd3be728f6488a2a78afb23a2051673a6bde97cf603f6292fdafea55cb2"
        ])
    );
}

#[test]
fn poseidon_bn254_t4_reference_kat() {
    let input: Vec<Bn254> = (1..=4).map(Bn254::from_u64).collect();
    let permutation = Poseidon::new(&POSEIDON_BN254_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input),
        expected_bn254(&[
            "0x13ff55f870c6a7f617515915d223e31e48282dfae157e36d6fa44a1b2a070306",
            "0x2def12f719945b6b58aa62a1f5ec9c46dec4f22b093be89f7a4894b85f3f2254",
            "0x249e4568c0f290f8c8960db6687b2c2c2b9b09867eed059b41a3d52c889ecc38",
            "0x2e177971715f689039d5923d841b0e1596b2fbfc5342cc67b8e5f4e68019d349",
        ])
    );
}
#[test]
fn poseidon2_bls12_t4_kat() {
    let permutation = Poseidon2::new(&POSEIDON2_BLS12_381_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x12bcaef1b22552236117b6bf614f681465d0531021698ae7ea2d42eb4f0ae928",
            "0xfd119eb1f08b8419afdb8f4e248fcb00391d2a71b0cdcb92192b6b376d4c6fe",
            "0x2f6dccc58c192115e6bbceaa16f967ae769cc75c84f457ff6e98940c7b22f6a6",
            "0x5279c9cfe15be66768ce00f03518ce1a6b9ed8c6fad7692a5b2e3475fb9419ab"
        ])
    );
}

#[test]
fn poseidon2_bn254_t4_gnark_crypto_kat() {
    let input: Vec<Bn254> = (1..=4).map(Bn254::from_u64).collect();
    let permutation = Poseidon2::new(&POSEIDON2_BN254_4_GNARK_CRYPTO_PARAMS);
    assert_eq!(
        permutation.permutation(&input),
        expected_bn254(&[
            "0x224785a48a72c75e2cbb698143e71d5d41bd89a2b9a7185871e39a54ce5785b1",
            "0x225bb800db22c4f4b09ace45cb484d42b0dd7dfe8708ee26aacde6f2c1fb2cb8",
            "0x1180f4260e60b4264c987b503075ea8374b53ed06c5145f8c21c2aadb5087d21",
            "0x16c877b5b9c04d873218804ccbf65d0eeb12db447f66c9ca26fec380055df7e9",
        ])
    );
}
#[test]
fn rescueprime_bls12_t4_kat() {
    let permutation = RescuePrime::new(&RESCUE_PRIME_BLS12_381_4_CURRENT_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x20a1320b43bcc6249a2d5ea73158c4d60142107f95566d74efffeea48041de39",
            "0x3228f6ff4c49533b0a07e825006c67cf256ade59f8846c9c0b6f7a0260e6c340",
            "0x16ca584a51a73a43e4d5df265c692a72ff5bc0b9a972f6f3a55d53c884dea534",
            "0x3a40ed2e2932c97fcb001dd745f423548712b4755d70f47d72eaa8ab267c5468"
        ])
    );
}
#[test]
fn anemoi_bls12_t4_kat() {
    let permutation = Anemoi::new(&ANEMOI_BLS12_381_4_CURRENT_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x3c09051f200539c064377225f55b112835bbaa6873a5aa74740586886472028f",
            "0x447c4e011098a61e3850beff200ac1190b7046d010c5ce5614778294748fddea",
            "0x10dd6f012115d2d105c35ea49824ef3908b467605516b8b40f069afe4ec9ccfa",
            "0x733d23ccc64ffa41ce416d32b13cbab1a9e2982e6ff460b8974a07b2d4e6e6cd"
        ])
    );
}
#[test]
fn griffin_bls12_t4_kat() {
    let permutation = Griffin::new(&GRIFFIN_BLS12_381_4_CURRENT_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x4ec9d82bc70461b509d6ea9e2ea4e2066a4f12e2d9dfab93a7294ed1dca6de59",
            "0xd813248c7b0ec4f7a7837b5c7a3bf483a48e29945bbff35bde888a997e46252",
            "0x3d3edd2a83287bc90ae4bd053e5bd4e0a75a3c8fa79abc3591f5729ed1629bce",
            "0x3a9be5012f515883dffb9f1517cec6b6316702055d78fe2224380f6baf63812e"
        ])
    );
}
#[test]
fn arion_bls12_t4_kat() {
    let permutation = Arion::new(&ARION_BLS12_381_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x6fc08546f8655be13d925451fc9fa08c7aedf7a6ce6dd28103384a60a740fa8c",
            "0x59ed9a0b1c71f64ec57a73a63e392d5f64c69a5968544fd964df48661bbc94d7",
            "0x605cdcdd514048a32f8324102eb5948c168c4aad04e998df6f6b8974fe0796f4",
            "0x5cdb55bec27e02252f1f0004c140040646afe73355f80b952fd38651a8bf8195"
        ])
    );
}
#[test]
fn polocolo_bls12_t4_kat() {
    let permutation = Polocolo::new(&POLOCOLO_BLS12_381_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x1069001b91b1ce7faf368634231d12fb7bd2b36d452e02ec15d0c812a75b6d13",
            "0x431778260687d961054ee973e656402fee36128ecd3efb7b6384a19a947d8beb",
            "0x491e1f769beba9f2f4d1d6360e5a68a385db205fa201b7dc40222cb946e645de",
            "0x1d03a3f1a0c0179e2230108e906f4fd8ff4ec085dfdf33733805f61b5c4b379e"
        ])
    );
}

#[test]
fn neptune_bls12_t4_kat() {
    let permutation = Neptune::new(&NEPTUNE_BLS12_381_4_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x16efa51cbd4b1ae3404e38551d53b132c5deda7eb899e0019ea3b986269742ec",
            "0x2958b61532cbdd2a0325d69e9d26bf16a50d7753a735dde5263aedfbfdf9b2fd",
            "0x6d9083403cac3312e2bfa9727a691ce756d13fd10ce675fc77f98511ba70e969",
            "0x34c90bdc1328ba7458573802927a7418eb4e568a90fc97b4a723e9cafc61603b"
        ])
    );
}

#[test]
fn skyscraper_bls12_t4_kat() {
    // Skyscraper has two extension-field lanes, so n=2 is t=4 over the base field.
    let permutation = Skyscraper::new(&SKYSCRAPER_BLS12_381_2_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x3a4bf00b4b4bddc463ebdd8e8fbc86c241b8ccf774a906df9b80d870b6dc55ef",
            "0x4403a64d0440af799ca0218fab4fb187ea87589af086cf20bf40cc7e1cc583fe",
            "0x321c05791d4e8e154b61685b9be51b0aff776011ca2d8e056251c2d4dd257a1",
            "0xe559dfab43557188bc923d6fc405e2d46ab54002280920211af8200a0a61c10"
        ])
    );
}

#[test]
fn grendel_bls12_t4_kat() {
    let permutation = Grendel::new(&GRENDEL_BLS12_381_4_CURRENT_PARAMS);
    assert_eq!(
        permutation.permutation(&input()),
        expected(&[
            "0x2e07b78265c9e344088d3b1a928c247ab59b4344c846b1bdb5e5f56cef1eb996",
            "0x2d7192cd99557716cf866e5ec35acfad4f398528ed7c96f1512ac372d0b526bc",
            "0x43b8112ee9d86e395f2ccc92b3362cfc5113da3f7bcd86eae52f745f9ab97689",
            "0x52d0162e77702bd0b41d5b6a2b7e8f3746af91f1488b11b74de8ad6ba23ab5d0"
        ])
    );
}

#[test]
fn bn254_exported_t4_profiles_smoke() {
    let input: Vec<Bn254> = (1..=4).map(Bn254::from_u64).collect();
    let outputs = [
        Gmimc::new(&crate::gmimc::instances::GMIMC_BN254_4_PARAMS).permutation(&input),
        Gmimc2::new(&crate::gmimc2::instances::GMIMC2_BN254_4_PARAMS).permutation(&input),
        Poseidon::new(&POSEIDON_BN254_4_PARAMS).permutation(&input),
        Poseidon2::new(&POSEIDON2_BN254_4_GNARK_CRYPTO_PARAMS).permutation(&input),
        Neptune::new(&crate::neptune::instances::NEPTUNE_BN254_4_PARAMS).permutation(&input),
        RescuePrime::new(&crate::rescueprime::instances::RESCUE_PRIME_BN254_4_CURRENT_PARAMS)
            .permutation(&input),
        Anemoi::new(&crate::anemoi::instances::ANEMOI_BN254_4_CURRENT_PARAMS).permutation(&input),
        Griffin::new(&crate::griffin::instances::GRIFFIN_BN254_4_CURRENT_PARAMS)
            .permutation(&input),
        Arion::new(&crate::arion::instances::ARION_BN254_4_PARAMS).permutation(&input),
        Skyscraper::new(&crate::skyscraper::instances::SKYSCRAPER_BN254_2_PARAMS)
            .permutation(&input),
        Grendel::new(&crate::grendel::instances::GRENDEL_BN254_4_CURRENT_PARAMS)
            .permutation(&input),
        Polocolo::new(&crate::polocolo::instances::POLOCOLO_BN254_4_PARAMS).permutation(&input),
    ];
    assert!(outputs.iter().all(|output| output.len() == 4));
}
