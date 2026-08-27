use core::cmp::Ordering;

pub trait MontyParams {
    const MODULUS: [u64; 4];
    const MU: u64;
    const R2: [u64; 4];
}

#[inline]
const fn carrying_add(lhs: u64, rhs: u64, carry: bool) -> (u64, bool) {
    let (a, c1) = lhs.overflowing_add(rhs);
    let (b, c2) = a.overflowing_add(carry as u64);
    (b, c1 | c2)
}

#[inline]
const fn borrowing_sub(lhs: u64, rhs: u64, borrow: bool) -> (u64, bool) {
    let (a, c1) = lhs.overflowing_sub(rhs);
    let (b, c2) = a.overflowing_sub(borrow as u64);
    (b, c1 | c2)
}

#[inline]
fn wrapping_add<const N: usize>(lhs: [u64; N], rhs: [u64; N]) -> ([u64; N], bool) {
    let mut carry = false;
    let mut output = [0u64; N];
    for i in 0..N {
        (output[i], carry) = carrying_add(lhs[i], rhs[i], carry);
    }
    (output, carry)
}

#[inline]
fn wrapping_sub<const N: usize>(lhs: [u64; N], rhs: [u64; N]) -> ([u64; N], bool) {
    let mut borrow = false;
    let mut output = [0u64; N];
    for i in 0..N {
        (output[i], borrow) = borrowing_sub(lhs[i], rhs[i], borrow);
    }
    (output, borrow)
}

#[inline]
fn cmp_limbs(lhs: [u64; 4], rhs: [u64; 4]) -> Ordering {
    for i in (0..4).rev() {
        if lhs[i] < rhs[i] {
            return Ordering::Less;
        }
        if lhs[i] > rhs[i] {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

#[inline]
fn mul_small(lhs: [u64; 4], rhs: u64) -> (u64, [u64; 4]) {
    let mut output = [0u64; 4];
    let mut acc = (lhs[0] as u128) * (rhs as u128);
    let out_0 = acc as u64;
    acc >>= 64;

    for i in 1..4 {
        acc += (lhs[i] as u128) * (rhs as u128);
        output[i - 1] = acc as u64;
        acc >>= 64;
    }
    output[3] = acc as u64;
    (out_0, output)
}

#[inline]
fn mul_small_and_acc(lhs: [u64; 4], rhs: u64, add: [u64; 4]) -> (u64, [u64; 4]) {
    let mut output = [0u64; 4];
    let mut acc = (lhs[0] as u128) * (rhs as u128) + (add[0] as u128);
    let out_0 = acc as u64;
    acc >>= 64;

    for i in 1..4 {
        acc += (lhs[i] as u128) * (rhs as u128) + (add[i] as u128);
        output[i - 1] = acc as u64;
        acc >>= 64;
    }
    output[3] = acc as u64;
    (out_0, output)
}

#[inline]
fn interleaved_monty_reduction<P: MontyParams>(acc0: u64, acc: [u64; 4]) -> [u64; 4] {
    let t = acc0.wrapping_mul(P::MU);
    let (_, u) = mul_small(P::MODULUS, t);

    let (sub, under) = wrapping_sub::<4>(acc, u);
    if under {
        let (sub_corr, _) = wrapping_add::<4>(sub, P::MODULUS);
        sub_corr
    } else {
        sub
    }
}

#[inline]
pub fn monty_mul<P: MontyParams>(lhs: [u64; 4], rhs: [u64; 4]) -> [u64; 4] {
    debug_assert!(cmp_limbs(lhs, P::MODULUS) == Ordering::Less);

    let (acc0, acc) = mul_small(lhs, rhs[0]);
    let res0 = interleaved_monty_reduction::<P>(acc0, acc);

    let (acc0, acc) = mul_small_and_acc(lhs, rhs[1], res0);
    let res1 = interleaved_monty_reduction::<P>(acc0, acc);

    let (acc0, acc) = mul_small_and_acc(lhs, rhs[2], res1);
    let res2 = interleaved_monty_reduction::<P>(acc0, acc);

    let (acc0, acc) = mul_small_and_acc(lhs, rhs[3], res2);
    interleaved_monty_reduction::<P>(acc0, acc)
}

#[inline]
pub fn add_mod<P: MontyParams>(lhs: [u64; 4], rhs: [u64; 4]) -> [u64; 4] {
    let (sum, carry) = wrapping_add::<4>(lhs, rhs);
    let (sum_minus, borrow) = wrapping_sub::<4>(sum, P::MODULUS);
    if carry || !borrow {
        sum_minus
    } else {
        sum
    }
}

#[inline]
pub fn sub_mod<P: MontyParams>(lhs: [u64; 4], rhs: [u64; 4]) -> [u64; 4] {
    let (diff, borrow) = wrapping_sub::<4>(lhs, rhs);
    if borrow {
        let (diff_corr, _) = wrapping_add::<4>(diff, P::MODULUS);
        diff_corr
    } else {
        diff
    }
}

#[inline]
pub fn to_monty<P: MontyParams>(value: [u64; 4]) -> [u64; 4] {
    monty_mul::<P>(value, P::R2)
}

pub fn reduce_raw<P: MontyParams>(mut value: [u64; 4]) -> [u64; 4] {
    while cmp_limbs(value, P::MODULUS) != Ordering::Less {
        value = sub_mod::<P>(value, P::MODULUS);
    }
    value
}

pub fn from_hex_to_limbs(s: &str) -> Option<[u64; 4]> {
    let s = s.trim();
    let s = s.strip_prefix("0x").unwrap_or(s);
    if s.is_empty() {
        return Some([0; 4]);
    }
    if s.len() > 64 {
        return None;
    }

    let padded_len = ((s.len() + 15) / 16) * 16;
    let mut padded = String::with_capacity(padded_len);
    for _ in 0..(padded_len - s.len()) {
        padded.push('0');
    }
    padded.push_str(s);

    let mut limbs = [0u64; 4];
    let mut idx = padded.len();
    let mut limb_idx = 0;
    while idx > 0 && limb_idx < 4 {
        let start = idx - 16;
        let chunk = &padded[start..idx];
        limbs[limb_idx] = u64::from_str_radix(chunk, 16).ok()?;
        idx = start;
        limb_idx += 1;
    }

    Some(limbs)
}
