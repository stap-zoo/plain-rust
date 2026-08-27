use core::marker::PhantomData;

pub trait Monty31Params {
    const PRIME: u32;
    const MONTY_BITS: u32 = 32;
    const MONTY_MASK: u32 = 0xFFFF_FFFF;
    const MONTY_MU: u32;
}

#[inline]
pub const fn to_monty<P: Monty31Params>(x: u32) -> u32 {
    (((x as u64) << P::MONTY_BITS) % (P::PRIME as u64)) as u32
}

#[inline]
pub const fn to_monty_64<P: Monty31Params>(x: u64) -> u32 {
    (((x as u128) << P::MONTY_BITS) % (P::PRIME as u128)) as u32
}

#[inline]
pub const fn from_monty<P: Monty31Params>(x: u32) -> u32 {
    monty_reduce::<P>(x as u64)
}

#[inline]
pub fn add<P: Monty31Params>(lhs: u32, rhs: u32) -> u32 {
    let mut sum = lhs + rhs;
    let (corr_sum, over) = sum.overflowing_sub(P::PRIME);
    if !over {
        sum = corr_sum;
    }
    sum
}

#[inline]
pub fn sub<P: Monty31Params>(lhs: u32, rhs: u32) -> u32 {
    let (mut diff, over) = lhs.overflowing_sub(rhs);
    let corr = if over { P::PRIME } else { 0 };
    diff = diff.wrapping_add(corr);
    diff
}

#[inline]
pub const fn monty_reduce<P: Monty31Params>(x: u64) -> u32 {
    let t = x.wrapping_mul(P::MONTY_MU as u64) & (P::MONTY_MASK as u64);
    let u = t * (P::PRIME as u64);

    let (x_sub_u, over) = x.overflowing_sub(u);
    let x_sub_u_hi = (x_sub_u >> P::MONTY_BITS) as u32;
    let corr = if over { P::PRIME } else { 0 };
    x_sub_u_hi.wrapping_add(corr)
}

#[derive(Clone, Copy, Default, Eq, PartialEq)]
#[repr(transparent)]
pub struct MontyField31<P: Monty31Params> {
    pub(crate) value: u32,
    _phantom: PhantomData<P>,
}

impl<P: Monty31Params> MontyField31<P> {
    #[inline(always)]
    pub const fn new(value: u32) -> Self {
        Self {
            value: to_monty::<P>(value),
            _phantom: PhantomData,
        }
    }

    #[inline(always)]
    pub const fn new_monty(value: u32) -> Self {
        Self {
            value,
            _phantom: PhantomData,
        }
    }

    #[inline(always)]
    pub const fn to_u32(self) -> u32 {
        from_monty::<P>(self.value)
    }

    #[inline]
    pub fn from_u64(value: u64) -> Self {
        Self {
            value: to_monty_64::<P>(value),
            _phantom: PhantomData,
        }
    }

    #[inline]
    pub fn add_assign(&mut self, rhs: &Self) {
        self.value = add::<P>(self.value, rhs.value);
    }

    #[inline]
    pub fn sub_assign(&mut self, rhs: &Self) {
        self.value = sub::<P>(self.value, rhs.value);
    }

    #[inline]
    pub fn mul_assign(&mut self, rhs: &Self) {
        self.value = monty_reduce::<P>((self.value as u64) * (rhs.value as u64));
    }
}
