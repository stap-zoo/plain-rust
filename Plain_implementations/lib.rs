#[path = "Anemoi/mod.rs"]
pub mod anemoi;
#[path = "Arion/mod.rs"]
pub mod arion;
pub mod fields;
#[path = "Gmimc/mod.rs"]
pub mod gmimc;
#[path = "Gmimc2/mod.rs"]
pub mod gmimc2;
#[path = "Grendel/mod.rs"]
pub mod grendel;
#[path = "Griffin/mod.rs"]
pub mod griffin;
#[path = "Monolith/mod.rs"]
pub mod monolith;
#[path = "Neptune/mod.rs"]
pub mod neptune;
pub mod plain_hashes;
#[path = "Polocolo/mod.rs"]
pub mod polocolo;
#[path = "Poseidon/mod.rs"]
pub mod poseidon;
#[path = "Poseidon2/mod.rs"]
pub mod poseidon2;
#[path = "Psquarehash/mod.rs"]
pub mod psquarehash;
#[path = "ReinforcedConcrete/mod.rs"]
pub mod reinforced_concrete;
#[path = "Rescueprime/mod.rs"]
pub mod rescueprime;
#[path = "Skyscraper/mod.rs"]
pub mod skyscraper;
#[path = "Tip4‘/mod.rs"]
pub mod tip4;
#[path = "Tip5/mod.rs"]
pub mod tip5;
mod utils;
#[path = "XHash/mod.rs"]
pub mod xhash;
// VisionMark32 removed: binary-field design, not in SoK benchmark scope.

#[cfg(test)]
mod paired_round_kats;
#[cfg(test)]
mod t4_kats;
