use super::monolith_params::{
    Monolith31Params, Monolith64Params, MonolithField32, MonolithField64,
};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct Monolith64<F: MonolithField64> {
    pub(crate) params: Arc<Monolith64Params<F>>,
}

#[derive(Clone, Debug)]
pub struct Monolith31<F: MonolithField32> {
    pub(crate) params: Arc<Monolith31Params<F>>,
}

// ====================================================================
// Monolith64 (Goldilocks: t=8, 12) — FFT-based circulant MDS from
// zkfriendlyhashzoo (mds_8.rs / mds_12.rs).  State is split into high
// and low 32-bit halves; each half is transformed independently, the
// two results are combined with a shift.
// ====================================================================

impl<F: MonolithField64> Monolith64<F> {
    pub fn new(params: &Arc<Monolith64Params<F>>) -> Self {
        Monolith64 {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut state = input.to_vec();
        self.concrete(&mut state, None);

        for rc in self.params.round_constants.iter() {
            self.bars(&mut state);
            self.bricks(&mut state);
            self.concrete(&mut state, Some(rc));
        }

        self.bars(&mut state);
        self.bricks(&mut state);
        self.concrete(&mut state, None);
        state
    }

    fn concrete(&self, state: &mut [F], rc: Option<&[F]>) {
        let t = state.len();
        let mut lo = [0u64; 12];
        let mut hi = [0u64; 12];
        for i in 0..t {
            let s = state[i].to_u64();
            hi[i] = s >> 32;
            lo[i] = (s as u32) as u64;
        }

        match t {
            8 => {
                let hi8: &[u64; 8] = hi[..8].try_into().unwrap();
                let lo8: &[u64; 8] = lo[..8].try_into().unwrap();
                let (l, h) = mds8_fft(lo8, hi8);
                for i in 0..8 {
                    lo[i] = l[i];
                    hi[i] = h[i];
                }
            }
            12 => {
                let hi12: &[u64; 12] = hi[..12].try_into().unwrap();
                let lo12: &[u64; 12] = lo[..12].try_into().unwrap();
                let (l, h) = mds12_fft(lo12, hi12);
                for i in 0..12 {
                    lo[i] = l[i];
                    hi[i] = h[i];
                }
            }
            _ => unreachable!(),
        }

        for i in 0..t {
            let mut out = F::from_u64(lo[i]);
            let mut hi_part = F::from_u64(hi[i]);
            hi_part.mul_assign(&F::from_u64(1u64 << 32));
            out.add_assign(&hi_part);
            if let Some(rc) = rc {
                out.add_assign(&rc[i]);
            }
            state[i] = out;
        }
    }

    fn bricks(&self, state: &mut [F]) {
        let prev = state.to_vec();
        for i in 1..state.len() {
            let mut sq = prev[i - 1].clone();
            sq.square();
            state[i].add_assign(&sq);
        }
    }

    fn bars(&self, state: &mut [F]) {
        for el in state.iter_mut().take(Monolith64Params::<F>::BARS) {
            let value = self.bar_u64_lookup(el.to_u64());
            *el = F::from_u64(value);
        }
    }

    fn bar_u64_lookup(&self, value: u64) -> u64 {
        let l1 = self.params.lookup[(value & 0xffff) as usize] as u64;
        let l2 = self.params.lookup[((value >> 16) & 0xffff) as usize] as u64;
        let l3 = self.params.lookup[((value >> 32) & 0xffff) as usize] as u64;
        let l4 = self.params.lookup[((value >> 48) & 0xffff) as usize] as u64;
        l1 | (l2 << 16) | (l3 << 32) | (l4 << 48)
    }
}

// ====================================================================
// Monolith31 (BabyBear/KoalaBear/Mersenne31: t=16, 24) — FFT-based
// circulant MDS for 32-bit fields.  Uses a Cooley-Tukey power-of-two
// decomposition from zkfriendlyhashzoo (mds_16.rs / mds_24.rs).
// ====================================================================

impl<F: MonolithField32> Monolith31<F> {
    pub fn new(params: &Arc<Monolith31Params<F>>) -> Self {
        Monolith31 {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut state = input.to_vec();
        self.concrete(&mut state, None);

        for rc in self.params.round_constants.iter() {
            self.bars(&mut state);
            self.bricks(&mut state);
            self.concrete(&mut state, Some(rc));
        }

        self.bars(&mut state);
        self.bricks(&mut state);
        self.concrete(&mut state, None);
        state
    }

    fn concrete(&self, state: &mut [F], rc: Option<&[F]>) {
        let t = state.len();
        let mut tmp = [0u64; 32];
        for i in 0..t {
            tmp[i] = state[i].to_u32() as u64;
        }

        match t {
            16 => {
                let r = fast_cyclomul16(*(&tmp[..16].try_into().unwrap()));
                for i in 0..16 {
                    let mut val = F::from_u64(r[i]);
                    if let Some(rc) = rc {
                        val.add_assign(&rc[i]);
                    }
                    state[i] = val;
                }
            }
            24 => {
                let mut buf = [0u64; 32];
                buf[0..24].copy_from_slice(&tmp[0..24]);
                let r = fast_cyclomul32(&buf);
                for i in 0..24 {
                    let mut val = F::from_u64(r[i]);
                    if let Some(rc) = rc {
                        val.add_assign(&rc[i]);
                    }
                    state[i] = val;
                }
            }
            _ => unreachable!(),
        }
    }

    fn bricks(&self, state: &mut [F]) {
        let prev = state.to_vec();
        for i in 1..state.len() {
            let mut sq = prev[i - 1].clone();
            sq.square();
            state[i].add_assign(&sq);
        }
    }

    fn bars(&self, state: &mut [F]) {
        for el in state.iter_mut().take(Monolith31Params::<F>::BARS) {
            let value = self.bar_u32_lookup(el.to_u32());
            *el = F::from_u64(value as u64);
        }
    }

    fn bar_u32_lookup(&self, value: u32) -> u32 {
        let low = self.params.lookup1[(value & 0xffff) as usize] as u32;
        let high = self.params.lookup2[(value >> 16) as usize] as u32;
        low | (high << 16)
    }
}

// ====================================================================
// Monolith64 FFT kernels — from zkfriendlyhashzoo mds_8.rs / mds_12.rs
// ====================================================================

/// t=8: MDS_FREQ constants for circulant first row [23,8,13,10,7,6,21,8]
const MDS8_BLOCK1: [i64; 2] = [16, 8];
const MDS8_BLOCK2: [(i64, i64); 2] = [(8, -4), (-1, 1)];
const MDS8_BLOCK3: [i64; 2] = [-1, 1];

/// t=12: MDS_FREQ constants for [7,23,8,26,13,10,9,7,6,22,21,8]
const MDS12_BLOCK1: [i64; 3] = [16, 8, 16];
const MDS12_BLOCK2: [(i64, i64); 3] = [(-1, 2), (-1, 1), (4, 8)];
const MDS12_BLOCK3: [i64; 3] = [-8, 1, 1];

#[inline(always)]
fn fft2_real(x: [u64; 2]) -> [i64; 2] {
    [(x[0] as i64 + x[1] as i64), (x[0] as i64 - x[1] as i64)]
}

#[inline(always)]
fn ifft2_real_unreduced(y: [i64; 2]) -> [u64; 2] {
    [(y[0] + y[1]) as u64, (y[0] - y[1]) as u64]
}

#[inline(always)]
fn fft4_real(x: [u64; 4]) -> (i64, (i64, i64), i64) {
    let [z0, z2] = fft2_real([x[0], x[2]]);
    let [z1, z3] = fft2_real([x[1], x[3]]);
    (z0 + z1, (z2, -z3), z0 - z1)
}

#[inline(always)]
fn ifft4_real_unreduced(y: (i64, (i64, i64), i64)) -> [u64; 4] {
    let z0 = y.0 + y.2;
    let z1 = y.0 - y.2;
    let z2 = y.1 .0;
    let z3 = -y.1 .1;
    let [x0, x2] = ifft2_real_unreduced([z0, z2]);
    let [x1, x3] = ifft2_real_unreduced([z1, z3]);
    [x0, x1, x2, x3]
}

// --- t=8 ---

fn mds8_fft(lo: &[u64; 8], hi: &[u64; 8]) -> ([u64; 8], [u64; 8]) {
    let lo = mds8_freq(*lo);
    let hi = mds8_freq(*hi);
    (lo, hi)
}

fn mds8_freq(state: [u64; 8]) -> [u64; 8] {
    let [s0, s1, s2, s3, s4, s5, s6, s7] = state;
    let (u0, u1, u2) = fft4_real([s0, s2, s4, s6]);
    let (u4, u5, u6) = fft4_real([s1, s3, s5, s7]);

    let [v0, v4] = mds8_block1([u0, u4], MDS8_BLOCK1);
    let [v1, v5] = mds8_block2([u1, u5], MDS8_BLOCK2);
    let [v2, v6] = mds8_block3([u2, u6], MDS8_BLOCK3);

    let [s0, s2, s4, s6] = ifft4_real_unreduced((v0, v1, v2));
    let [s1, s3, s5, s7] = ifft4_real_unreduced((v4, v5, v6));
    [s0, s1, s2, s3, s4, s5, s6, s7]
}

#[inline(always)]
fn mds8_block1(x: [i64; 2], y: [i64; 2]) -> [i64; 2] {
    [x[0] * y[0] + x[1] * y[1], x[0] * y[1] + x[1] * y[0]]
}

#[inline(always)]
fn mds8_block2(x: [(i64, i64); 2], y: [(i64, i64); 2]) -> [(i64, i64); 2] {
    let [(x0r, x0i), (x1r, x1i)] = x;
    let [(y0r, y0i), (y1r, y1i)] = y;
    let x0s = x0r + x0i;
    let x1s = x1r + x1i;
    let y0s = y0r + y0i;
    let y1s = y1r + y1i;
    let m0 = (x0r * y0r, x0i * y0i);
    let m1 = (x1r * y1r, x1i * y1i);
    let z0r = (m0.0 - m0.1) + (x1s * y1s - m1.0 - m1.1);
    let z0i = (x0s * y0s - m0.0 - m0.1) + (-m1.0 + m1.1);
    let m0 = (x0r * y1r, x0i * y1i);
    let m1 = (x1r * y0r, x1i * y0i);
    let z1r = (m0.0 - m0.1) + (m1.0 - m1.1);
    let z1i = (x0s * y1s - m0.0 - m0.1) + (x1s * y0s - m1.0 - m1.1);
    [(z0r, z0i), (z1r, z1i)]
}

#[inline(always)]
fn mds8_block3(x: [i64; 2], y: [i64; 2]) -> [i64; 2] {
    [x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]]
}

// --- t=12 ---

fn mds12_fft(lo: &[u64; 12], hi: &[u64; 12]) -> ([u64; 12], [u64; 12]) {
    let lo = mds12_freq(*lo);
    let hi = mds12_freq(*hi);
    (lo, hi)
}

fn mds12_freq(state: [u64; 12]) -> [u64; 12] {
    let [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11] = state;
    let (u0, u1, u2) = fft4_real([s0, s3, s6, s9]);
    let (u4, u5, u6) = fft4_real([s1, s4, s7, s10]);
    let (u8, u9, u10) = fft4_real([s2, s5, s8, s11]);

    let [v0, v4, v8] = mds12_block1([u0, u4, u8], MDS12_BLOCK1);
    let [v1, v5, v9] = mds12_block2([u1, u5, u9], MDS12_BLOCK2);
    let [v2, v6, v10] = mds12_block3([u2, u6, u10], MDS12_BLOCK3);

    let [s0, s3, s6, s9] = ifft4_real_unreduced((v0, v1, v2));
    let [s1, s4, s7, s10] = ifft4_real_unreduced((v4, v5, v6));
    let [s2, s5, s8, s11] = ifft4_real_unreduced((v8, v9, v10));
    [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11]
}

#[inline(always)]
fn mds12_block1(x: [i64; 3], y: [i64; 3]) -> [i64; 3] {
    [
        x[0] * y[0] + x[1] * y[2] + x[2] * y[1],
        x[0] * y[1] + x[1] * y[0] + x[2] * y[2],
        x[0] * y[2] + x[1] * y[1] + x[2] * y[0],
    ]
}

#[inline(always)]
fn mds12_block2(x: [(i64, i64); 3], y: [(i64, i64); 3]) -> [(i64, i64); 3] {
    let [(x0r, x0i), (x1r, x1i), (x2r, x2i)] = x;
    let [(y0r, y0i), (y1r, y1i), (y2r, y2i)] = y;
    let x0s = x0r + x0i;
    let x1s = x1r + x1i;
    let x2s = x2r + x2i;
    let y0s = y0r + y0i;
    let y1s = y1r + y1i;
    let y2s = y2r + y2i;

    let m0 = (x0r * y0r, x0i * y0i);
    let m1 = (x1r * y2r, x1i * y2i);
    let m2 = (x2r * y1r, x2i * y1i);
    let z0r = (m0.0 - m0.1) + (x1s * y2s - m1.0 - m1.1) + (x2s * y1s - m2.0 - m2.1);
    let z0i = (x0s * y0s - m0.0 - m0.1) + (-m1.0 + m1.1) + (-m2.0 + m2.1);
    let z0 = (z0r, z0i);

    let m0 = (x0r * y1r, x0i * y1i);
    let m1 = (x1r * y0r, x1i * y0i);
    let m2 = (x2r * y2r, x2i * y2i);
    let z1r = (m0.0 - m0.1) + (m1.0 - m1.1) + (x2s * y2s - m2.0 - m2.1);
    let z1i = (x0s * y1s - m0.0 - m0.1) + (x1s * y0s - m1.0 - m1.1) + (-m2.0 + m2.1);
    let z1 = (z1r, z1i);

    let m0 = (x0r * y2r, x0i * y2i);
    let m1 = (x1r * y1r, x1i * y1i);
    let m2 = (x2r * y0r, x2i * y0i);
    let z2r = (m0.0 - m0.1) + (m1.0 - m1.1) + (m2.0 - m2.1);
    let z2i = (x0s * y2s - m0.0 - m0.1) + (x1s * y1s - m1.0 - m1.1) + (x2s * y0s - m2.0 - m2.1);
    let z2 = (z2r, z2i);

    [z0, z1, z2]
}

#[inline(always)]
fn mds12_block3(x: [i64; 3], y: [i64; 3]) -> [i64; 3] {
    [
        x[0] * y[0] - x[1] * y[2] - x[2] * y[1],
        x[0] * y[1] + x[1] * y[0] - x[2] * y[2],
        x[0] * y[2] + x[1] * y[1] + x[2] * y[0],
    ]
}

// ====================================================================
// Monolith31 FFT kernels — from zkfriendlyhashzoo mds_16.rs / mds_24.rs
// Recursive Cooley-Tukey decomposition for power-of-two sizes.
// ====================================================================

// --- t=16: cyclomul16 → cyclomul8 → cyclomul4 → cyclomul2 ---

#[inline(always)]
fn fast_cyclomul16_31(f: &[u64; 16]) -> [u64; 16] {
    const N: usize = 8;
    let mut ff_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] as i64 + f[i + N] as i64;
        ff_hi[i] = f[i] as i64 - f[i + N] as i64;
    }
    let hh_lo = fast_cyclomul8_31(ff_lo);
    let hh_hi = complex_negacyclomul8_31(ff_hi);
    let mut hh = [0u64; 2 * N];
    for i in 0..N {
        hh[i] = ((hh_lo[i] + hh_hi[i]) >> 1) as u64;
        hh[i + N] = ((hh_lo[i] - hh_hi[i]) >> 1) as u64;
    }
    hh
}

fn fast_cyclomul16(f: &[u64; 16]) -> [u64; 16] {
    fast_cyclomul16_31(f)
}

#[inline(always)]
fn fast_cyclomul8_31(f: [i64; 8]) -> [i64; 8] {
    const N: usize = 4;
    let mut ff_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] + f[i + N];
        ff_hi[i] = f[i] - f[i + N];
    }
    let hh_lo = fast_cyclomul4_31(ff_lo);
    let hh_hi = complex_negacyclomul4_31(ff_hi);
    let mut hh = [0i64; 2 * N];
    for i in 0..N {
        hh[i] = (hh_lo[i] + hh_hi[i]) >> 1;
        hh[i + N] = (hh_lo[i] - hh_hi[i]) >> 1;
    }
    hh
}

#[inline(always)]
fn fast_cyclomul4_31(f: [i64; 4]) -> [i64; 4] {
    const N: usize = 2;
    let mut ff_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] + f[i + N];
        ff_hi[i] = f[i] - f[i + N];
    }
    let hh_lo = fast_cyclomul2_31(ff_lo);
    let hh_hi = complex_negacyclomul2_31(ff_hi);
    let mut hh = [0i64; 2 * N];
    for i in 0..N {
        hh[i] = (hh_lo[i] + hh_hi[i]) >> 1;
        hh[i + N] = (hh_lo[i] - hh_hi[i]) >> 1;
    }
    hh
}

#[inline(always)]
fn fast_cyclomul2_31(f: [i64; 2]) -> [i64; 2] {
    let ff_lo = f[0] + f[1];
    let ff_hi = f[0] - f[1];
    let hh_lo = ff_lo * 524757;
    let hh_hi = ff_hi * 52427;
    [(hh_lo + hh_hi) >> 1, (hh_lo - hh_hi) >> 1]
}

#[inline(always)]
fn complex_negacyclomul2_31(f: [i64; 2]) -> [i64; 2] {
    let f0 = (f[0], -f[1]);
    let g0 = (-12936, -26959);
    let h0 = (f0.0 * g0.0 - f0.1 * g0.1, f0.0 * g0.1 + f0.1 * g0.0);
    [h0.0, -h0.1]
}

#[inline(always)]
fn complex_negacyclomul4_31(f: [i64; 4]) -> [i64; 4] {
    let g0 = [(98878, 10562), (-74304, -44845)];
    let f0 = [(f[0], -f[2]), (f[1], -f[3])];
    let h0 = complex_karatsuba2_31(f0, g0);
    let mut h = [0i64; 7];
    for i in 0..3 {
        h[i] += h0[i].0;
        h[i + 2] -= h0[i].1;
    }
    let mut hh = [0i64; 4];
    for i in 0..4 {
        hh[i] += h[i];
    }
    for i in 4..7 {
        hh[i - 4] -= h[i];
    }
    hh
}

fn complex_karatsuba2_31(f: [(i64, i64); 2], g: [(i64, i64); 2]) -> [(i64, i64); 3] {
    let ff = (f[0].0 + f[1].0, f[0].1 + f[1].1);
    let gg = (g[0].0 + g[1].0, g[0].1 + g[1].1);
    let lo = (
        f[0].0 * g[0].0 - f[0].1 * g[0].1,
        f[0].0 * g[0].1 + f[0].1 * g[0].0,
    );
    let hi = (
        f[1].0 * g[1].0 - f[1].1 * g[1].1,
        f[1].0 * g[1].1 + f[1].1 * g[1].0,
    );
    let ffgg = (ff.0 * gg.0 - ff.1 * gg.1, ff.0 * gg.1 + ff.1 * gg.0);
    let li = (ffgg.0 - (lo.0 + hi.0), ffgg.1 - (lo.1 + hi.1));
    [lo, li, hi]
}

#[inline(always)]
fn complex_negacyclomul8_31(f: [i64; 8]) -> [i64; 8] {
    let mut f0 = [(0i64, 0i64); 4];
    let g0 = [
        (4451, 4567),
        (-26413, 16445),
        (-12601, -27067),
        (-7078, 5811),
    ];
    for i in 0..4 {
        f0[i] = (f[i], -f[4 + i]);
    }
    let h0 = complex_karatsuba4_31(f0, g0);
    let mut h = [0i64; 11];
    for i in 0..7 {
        h[i] += h0[i].0;
        h[i + 4] -= h0[i].1;
    }
    let mut hh = [0i64; 8];
    for i in 0..8 {
        hh[i] += h[i];
    }
    for i in 8..11 {
        hh[i - 8] -= h[i];
    }
    hh
}

fn complex_karatsuba4_31(f: [(i64, i64); 4], g: [(i64, i64); 4]) -> [(i64, i64); 7] {
    let ff = [
        (f[0].0 + f[2].0, f[0].1 + f[2].1),
        (f[1].0 + f[3].0, f[1].1 + f[3].1),
    ];
    let gg = [
        (g[0].0 + g[2].0, g[0].1 + g[2].1),
        (g[1].0 + g[3].0, g[1].1 + g[3].1),
    ];
    let lo = complex_karatsuba2_31([f[0], f[1]], [g[0], g[1]]);
    let hi = complex_karatsuba2_31([f[2], f[3]], [g[2], g[3]]);
    let mid = complex_karatsuba2_31(ff, gg);
    let mut li = [(0i64, 0i64); 3];
    for i in 0..3 {
        li[i].0 = mid[i].0 - (lo[i].0 + hi[i].0);
        li[i].1 = mid[i].1 - (lo[i].1 + hi[i].1);
    }
    let mut result = [(0i64, 0i64); 7];
    for i in 0..3 {
        result[i] = lo[i];
        result[2 + i].0 += li[i].0;
        result[2 + i].1 += li[i].1;
        result[4 + i] = hi[i];
    }
    result
}

// --- t=24: 32-point cyclomul (zero-padded to power of 2) ---
// 24-point uses 32-point: cyclomul32 → 16 → 8 → 4 → 2

#[inline(always)]
fn fast_cyclomul32(f: &[u64; 32]) -> [u64; 32] {
    const N: usize = 16;
    let mut ff_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] as i64 + f[i + N] as i64;
        ff_hi[i] = f[i] as i64 - f[i + N] as i64;
    }
    let hh_lo = fast_cyclomul16_24(ff_lo);
    let hh_hi = complex_negacyclomul16_24(ff_hi);
    let mut hh = [0u64; 2 * N];
    for i in 0..N {
        hh[i] = ((hh_lo[i] + hh_hi[i]) >> 1) as u64;
        hh[i + N] = ((hh_lo[i] - hh_hi[i]) >> 1) as u64;
    }
    hh
}

#[inline(always)]
fn fast_cyclomul16_24(f: [i64; 16]) -> [i64; 16] {
    const N: usize = 8;
    let mut ff_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] + f[i + N];
        ff_hi[i] = f[i] - f[i + N];
    }
    let hh_lo = fast_cyclomul8_24(ff_lo);
    let hh_hi = complex_negacyclomul8_24(ff_hi);
    let mut hh = [0i64; 2 * N];
    for i in 0..N {
        hh[i] = (hh_lo[i] + hh_hi[i]) >> 1;
        hh[i + N] = (hh_lo[i] - hh_hi[i]) >> 1;
    }
    hh
}

#[inline(always)]
fn fast_cyclomul8_24(f: [i64; 8]) -> [i64; 8] {
    const N: usize = 4;
    let mut ff_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] + f[i + N];
        ff_hi[i] = f[i] - f[i + N];
    }
    let hh_lo = fast_cyclomul4_24(ff_lo);
    let hh_hi = complex_negacyclomul4_24(ff_hi);
    let mut hh = [0i64; 2 * N];
    for i in 0..N {
        hh[i] = (hh_lo[i] + hh_hi[i]) >> 1;
        hh[i + N] = (hh_lo[i] - hh_hi[i]) >> 1;
    }
    hh
}

#[inline(always)]
fn fast_cyclomul4_24(f: [i64; 4]) -> [i64; 4] {
    const N: usize = 2;
    let mut ff_lo = [0i64; N];
    let mut ff_hi = [0i64; N];
    for i in 0..N {
        ff_lo[i] = f[i] + f[i + N];
        ff_hi[i] = f[i] - f[i + N];
    }
    let hh_lo = fast_cyclomul2_24(ff_lo);
    let hh_hi = complex_negacyclomul2_24(ff_hi);
    let mut hh = [0i64; 2 * N];
    for i in 0..N {
        hh[i] = (hh_lo[i] + hh_hi[i]) >> 1;
        hh[i + N] = (hh_lo[i] - hh_hi[i]) >> 1;
    }
    hh
}

#[inline(always)]
fn fast_cyclomul2_24(f: [i64; 2]) -> [i64; 2] {
    let ff_lo = f[0] + f[1];
    let ff_hi = f[0] - f[1];
    let hh_lo = ff_lo * 37460020068;
    let hh_hi = ff_hi * -2147483648;
    [(hh_lo + hh_hi) >> 1, (hh_lo - hh_hi) >> 1]
}

#[inline(always)]
fn complex_negacyclomul2_24(f: [i64; 2]) -> [i64; 2] {
    let f0 = (f[0], -f[1]);
    let g0 = (-32768, 4294934526);
    let h0 = (f0.0 * g0.0 - f0.1 * g0.1, f0.0 * g0.1 + f0.1 * g0.0);
    [h0.0, -h0.1]
}

#[inline(always)]
fn complex_negacyclomul4_24(f: [i64; 4]) -> [i64; 4] {
    let g0 = [(-1267653833, 1858422187), (-879829814, -289061460)];
    let f0 = [(f[0], -f[2]), (f[1], -f[3])];
    let h0 = complex_karatsuba2_24(f0, g0);
    let mut h = [0i64; 7];
    for i in 0..3 {
        h[i] += h0[i].0;
        h[i + 2] -= h0[i].1;
    }
    let mut hh = [0i64; 4];
    for i in 0..4 {
        hh[i] += h[i];
    }
    for i in 4..7 {
        hh[i - 4] -= h[i];
    }
    hh
}

fn complex_karatsuba2_24(f: [(i64, i64); 2], g: [(i64, i64); 2]) -> [(i64, i64); 3] {
    let ff = (f[0].0 + f[1].0, f[0].1 + f[1].1);
    let gg = (g[0].0 + g[1].0, g[0].1 + g[1].1);
    let lo = (
        f[0].0 * g[0].0 - f[0].1 * g[0].1,
        f[0].0 * g[0].1 + f[0].1 * g[0].0,
    );
    let hi = (
        f[1].0 * g[1].0 - f[1].1 * g[1].1,
        f[1].0 * g[1].1 + f[1].1 * g[1].0,
    );
    let ffgg = (ff.0 * gg.0 - ff.1 * gg.1, ff.0 * gg.1 + ff.1 * gg.0);
    let li = (ffgg.0 - (lo.0 + hi.0), ffgg.1 - (lo.1 + hi.1));
    [lo, li, hi]
}

#[inline(always)]
fn complex_negacyclomul8_24(f: [i64; 8]) -> [i64; 8] {
    let mut f0 = [(0i64, 0i64); 4];
    let g0 = [
        (-123384022, 726686671),
        (123384022, -1420796976),
        (-1316554499, -743987706),
        (743987706, -830929148),
    ];
    for i in 0..4 {
        f0[i] = (f[i], -f[4 + i]);
    }
    let h0 = complex_karatsuba4_24(f0, g0);
    let mut h = [0i64; 11];
    for i in 0..7 {
        h[i] += h0[i].0;
        h[i + 4] -= h0[i].1;
    }
    let mut hh = [0i64; 8];
    for i in 0..8 {
        hh[i] += h[i];
    }
    for i in 8..11 {
        hh[i - 8] -= h[i];
    }
    hh
}

fn complex_karatsuba4_24(f: [(i64, i64); 4], g: [(i64, i64); 4]) -> [(i64, i64); 7] {
    let ff = [
        (f[0].0 + f[2].0, f[0].1 + f[2].1),
        (f[1].0 + f[3].0, f[1].1 + f[3].1),
    ];
    let gg = [
        (g[0].0 + g[2].0, g[0].1 + g[2].1),
        (g[1].0 + g[3].0, g[1].1 + g[3].1),
    ];
    let lo = complex_karatsuba2_24([f[0], f[1]], [g[0], g[1]]);
    let hi = complex_karatsuba2_24([f[2], f[3]], [g[2], g[3]]);
    let mid = complex_karatsuba2_24(ff, gg);
    let mut li = [(0i64, 0i64); 3];
    for i in 0..3 {
        li[i].0 = mid[i].0 - (lo[i].0 + hi[i].0);
        li[i].1 = mid[i].1 - (lo[i].1 + hi[i].1);
    }
    let mut result = [(0i64, 0i64); 7];
    for i in 0..3 {
        result[i] = lo[i];
        result[2 + i].0 += li[i].0;
        result[2 + i].1 += li[i].1;
        result[4 + i] = hi[i];
    }
    result
}

#[inline(always)]
fn complex_negacyclomul16_24(f: [i64; 16]) -> [i64; 16] {
    let mut f0 = [(0i64, 0i64); 8];
    let g = [
        (-1653474029, -764368101),
        (-494009618, 1383115546),
        (-623486820, -787285939),
        (21828258, 107898302),
        (820828241, -1013121520),
        (1013121520, 1326655406),
        (-107898302, -21828258),
        (787285939, -1523996827),
    ];
    for i in 0..8 {
        f0[i] = (f[i], -f[8 + i]);
    }
    let h0 = complex_karatsuba8_24(f0, g);
    let mut h = [0i64; 23];
    for i in 0..15 {
        h[i] += h0[i].0;
        h[i + 8] -= h0[i].1;
    }
    let mut hh = [0i64; 16];
    for i in 0..16 {
        hh[i] += h[i];
    }
    for i in 16..23 {
        hh[i - 16] -= h[i];
    }
    hh
}

fn complex_karatsuba8_24(f: [(i64, i64); 8], g: [(i64, i64); 8]) -> [(i64, i64); 15] {
    let ff = [
        (f[0].0 + f[4].0, f[0].1 + f[4].1),
        (f[1].0 + f[5].0, f[1].1 + f[5].1),
        (f[2].0 + f[6].0, f[2].1 + f[6].1),
        (f[3].0 + f[7].0, f[3].1 + f[7].1),
    ];
    let gg = [
        (-832645788, -1777489621),
        (519111902, 2709770952),
        (-731385122, -809114197),
        (809114197, -1416098525),
    ];
    let lo = complex_karatsuba4_24([f[0], f[1], f[2], f[3]], [g[0], g[1], g[2], g[3]]);
    let hi = complex_karatsuba4_24([f[4], f[5], f[6], f[7]], [g[4], g[5], g[6], g[7]]);
    let mid = complex_karatsuba4_24(ff, gg);
    let mut li = [(0i64, 0i64); 7];
    for i in 0..7 {
        li[i].0 = mid[i].0 - (lo[i].0 + hi[i].0);
        li[i].1 = mid[i].1 - (lo[i].1 + hi[i].1);
    }
    let mut result = [(0i64, 0i64); 15];
    for i in 0..7 {
        result[i] = lo[i];
        result[4 + i].0 += li[i].0;
        result[4 + i].1 += li[i].1;
        result[8 + i] = hi[i];
    }
    result
}
