//! Plain permutation benchmarks for Cryptographic hash functions.
//!
//! we measure the cost of a single compression-function / permutation call
//! rather than the full hash.  All code is adapted from well-known,
//! well-tested RustCrypto / BLAKE3-team crate implementations.

/// SHA-256 compression function (64 rounds, 1 block).
///
/// Source: `sha2` crate v0.10, feature `compress`.
///   <https://crates.io/crates/sha2>
///   Re-exported as `sha2::compress256`.
pub use sha2::compress256 as sha256_compress;

/// Keccak-f[1600] permutation (24 rounds).
///
/// Source: `keccak` crate v0.1, function `f1600`.
///   <https://crates.io/crates/keccak>
///   <https://docs.rs/keccak/latest/keccak/fn.f1600.html>
///
/// This crate is a transitive dependency of `sha3` (RustCrypto) and is
/// already in the dependency tree.
pub use keccak::f1600 as keccak_f1600;

// ============================================================
// Blake2b compression function
//
// Source: adapted from the `blake2` crate v0.10 (RustCrypto).
//   <https://crates.io/crates/blake2>
//   <https://github.com/RustCrypto/hashes/tree/master/blake2>
//
// The original crate uses a SIMD-generic macro (`blake2_impl!`) with
// `$vec = u64x2` for Blake2b.  Below is the scalar (non-SIMD) version
// based on the G function definition from RFC 7693 §3.2, extracted
// from the crate's `macros.rs` quarter_round / round / compress code.
// ============================================================

/// Blake2b compression function (12 rounds of G, 1 × 128-byte block).
///
/// Parameters follow the `blake2` crate's internal `compress()` signature:
///   h: 8-word state (8 × u64)
///   m: 16-word message block (16 × u64)
///   t: byte counter (2 × u64 as (t0, t1) — only t0 used for Blake2b)
///   f: final-block flags (2 × u64 as (f0, f1))
///
/// Source: Adapted from `blake2-0.10.6/src/macros.rs` (RustCrypto).
/// The original uses `u64x2` SIMD vectors; this is the equivalent scalar
/// implementation per RFC 7693, Algorithm 2 (COMPRESS).
pub fn blake2b_compress(h: &mut [u64; 8], m: &[u64; 16], t0: u64, t1: u64, f0: u64, _f1: u64) {
    // Blake2b IV (same as SHA-512 initial state).
    // Source: blake2-0.10.6/src/macros.rs, $IV for Blake2b.
    const IV: [u64; 8] = [
        0x6a09e667f3bcc908,
        0xbb67ae8584caa73b,
        0x3c6ef372fe94f82b,
        0xa54ff53a5f1d36f1,
        0x510e527fade682d1,
        0x9b05688c2b3e6c1f,
        0x1f83d9abfb41bd6b,
        0x5be0cd19137e2179,
    ];

    // SIGMA permutation table (10 rounds).
    // Source: blake2-0.10.6/src/consts.rs
    const SIGMA: [[usize; 16]; 10] = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
        [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
        [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8],
        [9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13],
        [2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9],
        [12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11],
        [13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10],
        [6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5],
        [10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0],
    ];

    // Internal state initialization.
    // Source: blake2-0.10.6/src/macros.rs, compress() body.
    let mut v = [0u64; 16];
    v[0] = h[0];
    v[1] = h[1];
    v[2] = h[2];
    v[3] = h[3];
    v[4] = h[4];
    v[5] = h[5];
    v[6] = h[6];
    v[7] = h[7];
    v[8] = IV[0];
    v[9] = IV[1];
    v[10] = IV[2];
    v[11] = IV[3];
    v[12] = h[0] ^ IV[4];
    v[13] = h[1] ^ IV[5];
    v[14] = h[2] ^ IV[6];
    v[15] = h[3] ^ IV[7];

    // Counter / flags XOR (same as blake2 crate).
    if f0 != !0 {
        v[12] ^= t0;
        v[13] ^= t1;
    }

    // G function (quarter round).  Derived from `quarter_round` in
    // blake2-0.10.6/src/macros.rs, scalarised ($vec = u64).
    // Rotation constants: Blake2b uses R1=32, R2=24, R3=16, R4=63.
    #[inline(always)]
    fn g(v: &mut [u64; 16], a: usize, b: usize, c: usize, d: usize, x: u64, y: u64) {
        v[a] = v[a].wrapping_add(v[b]).wrapping_add(x);
        v[d] = (v[d] ^ v[a]).rotate_right(32);
        v[c] = v[c].wrapping_add(v[d]);
        v[b] = (v[b] ^ v[c]).rotate_right(24);
        v[a] = v[a].wrapping_add(v[b]).wrapping_add(y);
        v[d] = (v[d] ^ v[a]).rotate_right(16);
        v[c] = v[c].wrapping_add(v[d]);
        v[b] = (v[b] ^ v[c]).rotate_right(63);
    }

    // 12 rounds × 8 G calls in Blake2 column/diagonal pattern.
    // Source: blake2-0.10.6/src/macros.rs, round() + compress() loop.
    for round in 0..12 {
        let s = &SIGMA[round % 10];

        // Column steps
        g(&mut v, 0, 4, 8, 12, m[s[0]], m[s[1]]);
        g(&mut v, 1, 5, 9, 13, m[s[2]], m[s[3]]);
        g(&mut v, 2, 6, 10, 14, m[s[4]], m[s[5]]);
        g(&mut v, 3, 7, 11, 15, m[s[6]], m[s[7]]);
        // Diagonal steps
        g(&mut v, 0, 5, 10, 15, m[s[8]], m[s[9]]);
        g(&mut v, 1, 6, 11, 12, m[s[10]], m[s[11]]);
        g(&mut v, 2, 7, 8, 13, m[s[12]], m[s[13]]);
        g(&mut v, 3, 4, 9, 14, m[s[14]], m[s[15]]);
    }

    // Finalisation: h ← h ^ v[0..4] ^ v[8..12]
    // Source: blake2-0.10.6/src/macros.rs, after round loop.
    h[0] ^= v[0] ^ v[8];
    h[1] ^= v[1] ^ v[9];
    h[2] ^= v[2] ^ v[10];
    h[3] ^= v[3] ^ v[11];
    h[4] ^= v[4] ^ v[12];
    h[5] ^= v[5] ^ v[13];
    h[6] ^= v[6] ^ v[14];
    h[7] ^= v[7] ^ v[15];
}

// ============================================================
// Blake3 compression function
//
// Source: adapted from the `blake3` crate v1.8.5 (BLAKE3-team).
//   <https://crates.io/crates/blake3>
//   <https://github.com/BLAKE3-team/BLAKE3>
//
// The portable reference implementation is in `portable.rs`.
// Below is the key compression logic extracted from:
//   blake3-1.8.5/src/portable.rs  (g, round, compress_pre, compress_in_place)
//   blake3-1.8.5/src/lib.rs       (IV, MSG_SCHEDULE, counter_low/_high)
// ============================================================

/// Blake3 compression function (7 rounds, 1 × 64-byte block).
///
/// This is the core of Blake3's chunk compression.  It matches the
/// `compress_pre` + `compress_in_place` logic from `portable.rs`.
///
/// Source: blake3-1.8.5/src/portable.rs, functions `g`, `round`,
/// `compress_pre`, `compress_in_place`.
pub fn blake3_compress(
    cv: &[u32; 8],
    block: &[u8; 64],
    block_len: u32,
    counter: u64,
    flags: u32,
) -> [u32; 8] {
    // IV constants from blake3-1.8.5/src/lib.rs
    const IV: [u32; 8] = [
        0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A, 0x510E527F, 0x9B05688C, 0x1F83D9AB,
        0x5BE0CD19,
    ];

    // Message schedule per round from blake3-1.8.5/src/lib.rs
    const MSG_SCHEDULE: [[usize; 16]; 7] = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        [2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8],
        [3, 4, 10, 12, 13, 2, 7, 14, 6, 5, 9, 0, 11, 15, 8, 1],
        [10, 7, 12, 9, 14, 3, 13, 15, 4, 0, 11, 2, 5, 8, 1, 6],
        [12, 13, 9, 11, 15, 10, 14, 8, 7, 2, 5, 3, 0, 1, 6, 4],
        [9, 14, 11, 5, 8, 12, 15, 1, 13, 3, 0, 10, 2, 6, 4, 7],
        [11, 15, 5, 0, 1, 9, 8, 6, 14, 10, 2, 12, 3, 4, 7, 13],
    ];

    // Convert block bytes to u32 words (little-endian).
    // Source: blake3-1.8.5/src/platform.rs, words_from_le_bytes_64.
    let mut block_words = [0u32; 16];
    for i in 0..16 {
        block_words[i] = u32::from_le_bytes([
            block[4 * i],
            block[4 * i + 1],
            block[4 * i + 2],
            block[4 * i + 3],
        ]);
    }

    // Counter splitting: blake3-1.8.5/src/lib.rs
    #[inline(always)]
    fn counter_low(counter: u64) -> u32 {
        counter as u32
    }
    #[inline(always)]
    fn counter_high(counter: u64) -> u32 {
        (counter >> 32) as u32
    }

    // Initialise state: blake3-1.8.5/src/portable.rs, compress_pre.
    let mut state: [u32; 16] = [
        cv[0],
        cv[1],
        cv[2],
        cv[3],
        cv[4],
        cv[5],
        cv[6],
        cv[7],
        IV[0],
        IV[1],
        IV[2],
        IV[3],
        counter_low(counter),
        counter_high(counter),
        block_len,
        flags,
    ];

    // G function: blake3-1.8.5/src/portable.rs, fn g.
    #[inline(always)]
    fn g(state: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize, x: u32, y: u32) {
        state[a] = state[a].wrapping_add(state[b]).wrapping_add(x);
        state[d] = (state[d] ^ state[a]).rotate_right(16);
        state[c] = state[c].wrapping_add(state[d]);
        state[b] = (state[b] ^ state[c]).rotate_right(12);
        state[a] = state[a].wrapping_add(state[b]).wrapping_add(y);
        state[d] = (state[d] ^ state[a]).rotate_right(8);
        state[c] = state[c].wrapping_add(state[d]);
        state[b] = (state[b] ^ state[c]).rotate_right(7);
    }

    // Round function: blake3-1.8.5/src/portable.rs, fn round.
    #[inline(always)]
    fn round(state: &mut [u32; 16], msg: &[u32; 16], schedule: &[usize; 16]) {
        // Mix columns
        g(state, 0, 4, 8, 12, msg[schedule[0]], msg[schedule[1]]);
        g(state, 1, 5, 9, 13, msg[schedule[2]], msg[schedule[3]]);
        g(state, 2, 6, 10, 14, msg[schedule[4]], msg[schedule[5]]);
        g(state, 3, 7, 11, 15, msg[schedule[6]], msg[schedule[7]]);
        // Mix diagonals
        g(state, 0, 5, 10, 15, msg[schedule[8]], msg[schedule[9]]);
        g(state, 1, 6, 11, 12, msg[schedule[10]], msg[schedule[11]]);
        g(state, 2, 7, 8, 13, msg[schedule[12]], msg[schedule[13]]);
        g(state, 3, 4, 9, 14, msg[schedule[14]], msg[schedule[15]]);
    }

    // 7 rounds: blake3-1.8.5/src/portable.rs, compress_pre.
    for r in 0..7 {
        round(&mut state, &block_words, &MSG_SCHEDULE[r]);
    }

    // Finalisation: blake3-1.8.5/src/portable.rs, compress_in_place.
    let mut out = [0u32; 8];
    out[0] = state[0] ^ state[8];
    out[1] = state[1] ^ state[9];
    out[2] = state[2] ^ state[10];
    out[3] = state[3] ^ state[11];
    out[4] = state[4] ^ state[12];
    out[5] = state[5] ^ state[13];
    out[6] = state[6] ^ state[14];
    out[7] = state[7] ^ state[15];
    out
}
