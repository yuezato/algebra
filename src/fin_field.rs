use crate::field::*;
use crate::univariate_polynomial::*;
use std::convert::TryInto;
use std::convert::{From, Into};
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait FiniteField: Field {
    // 全要素を列挙する
    fn enumerate() -> Vec<Self>;

    // 体の要素として埋め込めるバイト数
    const BYTE_SIZE: usize;
    const CARDINALITY: usize;
    type RepType;

    fn from_bytes(v: &[u8]) -> Self;

    fn to_byte(&self, idx: usize) -> u8;

    fn to_usize(&self) -> usize;

    fn get(v: &[u8], idx: usize) -> Self {
        debug_assert!(v.len() % Self::BYTE_SIZE == 0);
        debug_assert!(v.len() >= (idx + 1) * Self::BYTE_SIZE);

        Self::from_bytes(&v[idx * Self::BYTE_SIZE..(idx + 1) * Self::BYTE_SIZE])
    }

    fn put(&self, v: &mut [u8], idx: usize) {
        debug_assert!(v.len() % Self::BYTE_SIZE == 0);
        debug_assert!(v.len() >= (idx + 1) * Self::BYTE_SIZE);

        for i in 0..Self::BYTE_SIZE {
            v[idx * Self::BYTE_SIZE + i] = self.to_byte(i);
        }
    }

    // v1 += k*v2;
    fn mul_then_add(k: Self, v1: &mut [u8], v2: &[u8]);

    // v1 += k*w2;
    fn mul_then_add2(mul_table: &[Self], v1: &mut [u8], v2: &[u8]);
}

pub trait HasPrimitiveElement: Copy {
    const PRIMITIVE_ELEMENT: Self;
}

/*
 * The section of GF(2)
 */

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(non_camel_case_types)]
pub struct GF_2 {
    // true = 1, false = 0
    value: bool,
}

/*
 *  + 0 1   * 0 1
 *  0 0 1   0 0 0
 *  1 1 0   1 0 1
 *
 *  + == bit XOR
 *  * == bit AND
 */
impl GF_2 {
    const ZERO: GF_2 = GF_2 { value: false };
    const ONE: GF_2 = GF_2 { value: true };

    fn mul_inv(&self) -> GF_2 {
        GF_2::ONE
    }
}

impl Add for GF_2 {
    type Output = GF_2;

    // XOR
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, rhs: GF_2) -> GF_2 {
        let value = self.value ^ rhs.value;
        GF_2 { value }
    }
}

impl Neg for GF_2 {
    type Output = GF_2;

    fn neg(self) -> GF_2 {
        self
    }
}

impl Sub for GF_2 {
    type Output = GF_2;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: GF_2) -> GF_2 {
        self + (-rhs)
    }
}

impl Mul for GF_2 {
    type Output = GF_2;

    // AND
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: GF_2) -> GF_2 {
        let value = self.value & rhs.value;
        GF_2 { value }
    }
}

impl Div for GF_2 {
    type Output = GF_2;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: GF_2) -> GF_2 {
        self * rhs.mul_inv()
    }
}

impl Field for GF_2 {
    const ZERO: GF_2 = GF_2::ZERO;
    const ONE: GF_2 = GF_2::ONE;

    fn mul_inv(&self) -> Self {
        self.mul_inv()
    }
}

impl FiniteField for GF_2 {
    fn enumerate() -> Vec<Self> {
        vec![GF_2::ZERO, GF_2::ONE]
    }

    const BYTE_SIZE: usize = 0;
    const CARDINALITY: usize = 2;
    type RepType = bool;

    fn from_bytes(v: &[u8]) -> GF_2 {
        debug_assert!(v.len() == 1);

        if v[0] == 1 {
            GF_2::ONE
        } else if v[0] == 0 {
            GF_2::ZERO
        } else {
            panic!("invalid value");
        }
    }

    fn to_byte(&self, idx: usize) -> u8 {
        debug_assert!(idx == 0);

        if !self.value {
            0
        } else {
            1
        }
    }

    fn to_usize(&self) -> usize {
        if !self.value {
            0
        } else {
            1
        }
    }

    fn mul_then_add(_k: Self, _v1: &mut [u8], _v2: &[u8]) {
        panic!("undefined");
    }

    fn mul_then_add2(_mul_table: &[Self], _v1: &mut [u8], _v2: &[u8]) {
        panic!("undefined");
    }
}

/*
 * The section of GF(2^8)
 */

lazy_static! {
    // use x^8 + x^4 + x^3 + x^2 + 1.
    pub static ref GF_2_8_IMPL: GF_2_8_impl = GF_2_8_impl::new(
        Poly::from_vec( vec![
            (8, GF_2::ONE), (4, GF_2::ONE), (3, GF_2::ONE), (2, GF_2::ONE), (0, GF_2::ONE)
        ])
    );
}

#[allow(non_camel_case_types)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct GF_2_8(u8);

impl GF_2_8 {
    // 1 \alpha^1 + 0 \alpha^0
    pub const PRIMITIVE_ROOT: GF_2_8 = GF_2_8(0b10);
}

impl HasPrimitiveElement for GF_2_8 {
    const PRIMITIVE_ELEMENT: GF_2_8 = GF_2_8(0b10);
}

impl Field for GF_2_8 {
    // 0
    const ZERO: GF_2_8 = GF_2_8(0);

    // \alpha^0
    const ONE: GF_2_8 = GF_2_8(1);

    fn mul_inv(&self) -> GF_2_8 {
        GF_2_8_IMPL.mul_inv(*self)
    }
}

impl FiniteField for GF_2_8 {
    fn enumerate() -> Vec<Self> {
        (0u8..=0xff).map(|v| v.into()).collect()
    }

    const BYTE_SIZE: usize = 1;
    const CARDINALITY: usize = 0x100;
    type RepType = u8;

    fn from_bytes(v: &[u8]) -> Self {
        debug_assert!(v.len() == 1);

        GF_2_8(v[0])
    }

    fn to_byte(&self, idx: usize) -> u8 {
        debug_assert!(idx == 0);

        self.0
    }

    fn to_usize(&self) -> usize {
        self.0.into()
    }

    fn mul_then_add(k: Self, v1: &mut [u8], v2: &[u8]) {
        for i in 0..v1.len() {
            let x: u8 = (k * v2[i].into()).into();
            v1[i] ^= x;
        }
    }

    fn mul_then_add2(mul_table: &[Self], v1: &mut [u8], v2: &[u8]) {
        for i in 0..v1.len() {
            let x = mul_table[v2[i] as usize];
            let x: u8 = x.into();
            v1[i] ^= x;
        }
    }
}

impl From<u8> for GF_2_8 {
    fn from(v: u8) -> Self {
        GF_2_8(v)
    }
}

impl From<GF_2_8> for u8 {
    fn from(v: GF_2_8) -> Self {
        v.0
    }
}

impl From<Poly<GF_2>> for GF_2_8 {
    fn from(p: Poly<GF_2>) -> Self {
        let mut v: u8 = 0;

        for (deg, coef) in p.iter() {
            if *coef == GF_2::ONE {
                v |= 1 << (*deg)
            }
        }

        v.into()
    }
}

impl Mul<GF_2_8> for GF_2_8 {
    type Output = GF_2_8;

    fn mul(self, rhs: GF_2_8) -> GF_2_8 {
        GF_2_8_IMPL.mul(self, rhs)
    }
}

impl Add<GF_2_8> for GF_2_8 {
    type Output = GF_2_8;

    fn add(self, rhs: GF_2_8) -> GF_2_8 {
        GF_2_8_IMPL.add(self, rhs)
    }
}

impl Neg for GF_2_8 {
    type Output = GF_2_8;

    fn neg(self) -> GF_2_8 {
        GF_2_8_IMPL.add_inv(self)
    }
}

impl Sub<GF_2_8> for GF_2_8 {
    type Output = GF_2_8;

    fn sub(self, rhs: GF_2_8) -> GF_2_8 {
        self + (-rhs)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div<GF_2_8> for GF_2_8 {
    type Output = GF_2_8;

    fn div(self, rhs: GF_2_8) -> GF_2_8 {
        self * GF_2_8_IMPL.mul_inv(rhs)
    }
}

#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
pub struct GF_2_8_impl {
    // primitive polynomial
    ppoly: Poly<GF_2>,
    psi: Vec<GF_2_8>,
    phi: Vec<u8>,
}

impl GF_2_8_impl {
    const MAX_EXP: u8 = 0xff - 1;
    const MODULO: u8 = 0xff;
    const ORDER: u16 = 0x100;

    pub fn zero() -> GF_2_8 {
        0.into()
    }

    pub fn one() -> GF_2_8 {
        1.into()
    }

    pub fn new(ppoly: Poly<GF_2>) -> GF_2_8_impl {
        debug_assert!(ppoly.degree() == Some(8));

        let mut psi: Vec<GF_2_8> = vec![0.into(); (GF_2_8_impl::MAX_EXP + 1) as usize];
        let mut phi: Vec<u8> = vec![0; GF_2_8_impl::ORDER as usize];
        let mut p = None;

        for i in 0u8..=GF_2_8_impl::MAX_EXP {
            if let Some(p_) = p {
                p = Some(p_ * Poly::<GF_2>::from_mono(1, GF_2::ONE));
            } else {
                p = Some(Poly::<GF_2>::from_mono(0, GF_2::ONE))
            }

            let reduced = &(p.unwrap()) % &ppoly;
            p = Some(reduced.clone());

            let rep: GF_2_8 = reduced.into();
            let bin_rep: u8 = rep.into();

            psi[i as usize] = rep;
            phi[bin_rep as usize] = i;
        }

        GF_2_8_impl { ppoly, psi, phi }
    }

    pub fn ppoly(&self) -> &Poly<GF_2> {
        &self.ppoly
    }

    pub fn mul(&self, p: GF_2_8, q: GF_2_8) -> GF_2_8 {
        let p = u8::from(p);
        let q = u8::from(q);

        if p == 0 || q == 0 {
            return 0.into();
        }

        let i: u16 = self.phi[p as usize].into();
        let j: u16 = self.phi[q as usize].into();

        let i_j = (i + j) % (GF_2_8_impl::MODULO as u16);
        self.psi[i_j as usize]
    }

    pub fn mul_inv(&self, p: GF_2_8) -> GF_2_8 {
        let p = u8::from(p);

        debug_assert!(p != 0);

        let i = self.phi[p as usize];
        let inv = (GF_2_8_impl::MODULO - i) % GF_2_8_impl::MODULO;
        self.psi[inv as usize]
    }

    pub fn add(&self, p: GF_2_8, q: GF_2_8) -> GF_2_8 {
        (u8::from(p) ^ u8::from(q)).into()
    }

    pub fn add_inv(&self, p: GF_2_8) -> GF_2_8 {
        p
    }
}

/*
 * The section of GF(2^16)
 */

lazy_static! {
    // use x^16 + x^12 + x^3 + x^1 + 1.
    pub static ref GF_2_16_IMPL: GF_2_16 = GF_2_16::new(
        Poly::from_vec( vec![
            (16, GF_2::ONE), (12, GF_2::ONE), (3, GF_2::ONE), (1, GF_2::ONE), (0, GF_2::ONE)
        ])
    );
}

#[allow(non_camel_case_types)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct GF_2_16_Val(u16);

impl GF_2_16_Val {
    // 1 \alpha^1 + 0 \alpha^0
    pub const PRIMITIVE_ROOT: GF_2_16_Val = GF_2_16_Val(0b10);
}

impl HasPrimitiveElement for GF_2_16_Val {
    const PRIMITIVE_ELEMENT: GF_2_16_Val = GF_2_16_Val(0b10);
}

impl Field for GF_2_16_Val {
    // 0
    const ZERO: GF_2_16_Val = GF_2_16_Val(0);

    // \alpha^0
    const ONE: GF_2_16_Val = GF_2_16_Val(1);

    fn mul_inv(&self) -> GF_2_16_Val {
        GF_2_16_IMPL.mul_inv(*self)
    }
}

impl FiniteField for GF_2_16_Val {
    fn enumerate() -> Vec<Self> {
        (0u16..=0xffff).map(|v| v.into()).collect()
    }

    const BYTE_SIZE: usize = 2;
    const CARDINALITY: usize = 0x1_00_00;
    type RepType = u16;

    fn from_bytes(v: &[u8]) -> Self {
        debug_assert!(v.len() == 2);

        let r: u16 = (u16::from(v[1]) << 8) | u16::from(v[0]);
        GF_2_16_Val(r)
    }

    fn to_byte(&self, idx: usize) -> u8 {
        debug_assert!(idx == 0 || idx == 1);

        let v: u8 = ((self.0 >> (8 * idx)) & 0xff).try_into().unwrap();
        v
    }

    fn to_usize(&self) -> usize {
        self.0.into()
    }

    fn mul_then_add(k: Self, v1: &mut [u8], v2: &[u8]) {
        for i in (0..v1.len()).step_by(Self::BYTE_SIZE) {
            let x = k * Self::from_bytes(&v2[i..i + 2]);
            v1[i] ^= x.to_byte(0);
            v1[i + 1] ^= x.to_byte(1);
        }
    }

    fn mul_then_add2(mul_table: &[Self], v1: &mut [u8], v2: &[u8]) {
        for i in (0..v1.len()).step_by(Self::BYTE_SIZE) {
            let idx: u16 = Self::from_bytes(&v2[i..i + 2]).into();
            let x: Self = mul_table[idx as usize];
            v1[i] ^= x.to_byte(0);
            v1[i + 1] ^= x.to_byte(1);
        }
    }
}

impl From<u16> for GF_2_16_Val {
    fn from(v: u16) -> Self {
        GF_2_16_Val(v)
    }
}

impl From<GF_2_16_Val> for u16 {
    fn from(v: GF_2_16_Val) -> Self {
        v.0
    }
}

impl From<Poly<GF_2>> for GF_2_16_Val {
    fn from(p: Poly<GF_2>) -> Self {
        let mut v: u16 = 0;

        for (deg, coef) in p.iter() {
            if *coef == GF_2::ONE {
                v |= 1 << (*deg)
            }
        }

        v.into()
    }
}

impl Mul<GF_2_16_Val> for GF_2_16_Val {
    type Output = GF_2_16_Val;

    fn mul(self, rhs: GF_2_16_Val) -> GF_2_16_Val {
        GF_2_16_IMPL.mul(self, rhs)
    }
}

impl Add<GF_2_16_Val> for GF_2_16_Val {
    type Output = GF_2_16_Val;

    fn add(self, rhs: GF_2_16_Val) -> GF_2_16_Val {
        GF_2_16_IMPL.add(self, rhs)
    }
}

impl Neg for GF_2_16_Val {
    type Output = GF_2_16_Val;

    fn neg(self) -> GF_2_16_Val {
        GF_2_16_IMPL.add_inv(self)
    }
}

impl Sub<GF_2_16_Val> for GF_2_16_Val {
    type Output = GF_2_16_Val;

    fn sub(self, rhs: GF_2_16_Val) -> GF_2_16_Val {
        self + (-rhs)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div<GF_2_16_Val> for GF_2_16_Val {
    type Output = GF_2_16_Val;

    fn div(self, rhs: GF_2_16_Val) -> GF_2_16_Val {
        self * GF_2_16_IMPL.mul_inv(rhs)
    }
}

#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
pub struct GF_2_16 {
    // primitive polynomial
    ppoly: Poly<GF_2>,

    /// psi: \alpha^i -> GF_2_16_Val
    /// psi: Map<u16, GF_2_16_Val>
    psi: Vec<GF_2_16_Val>,

    /// phi: GF_2_16_Val -> 0 or \alpha^i
    /// phi: Map<GF_2_16_Val, u16>
    phi: Vec<u16>,
}

impl GF_2_16 {
    /*
     * (\alpha^0, \alpha^1, ..., \alpha^max_exp) is
    `* the multiplicative roup of GF(2^16).
     * Since \alpha is a primitive root, \alpha^{\max_exp+1} = 1.
     */
    const MAX_EXP: u16 = 0xffff - 1;

    /// Since modulo = max_exp + 1,
    /// for any p % `modulo` = 0, \alpha^p = 1.
    const MODULO: u16 = 0xffff;

    // The order of GF(2^16).
    const ORDER: u32 = 0x10000;

    pub fn zero() -> GF_2_16_Val {
        0.into()
    }

    pub fn one() -> GF_2_16_Val {
        1.into()
    }

    pub fn new(ppoly: Poly<GF_2>) -> GF_2_16 {
        assert!(ppoly.degree() == Some(16));

        // psi : \alpha^i -> reduced poly
        // Therefore, psi[0] is the (1 % P).
        // psi[0] is not the 0.
        let mut psi: Vec<GF_2_16_Val> =
            vec![0.into(); GF_2_16::MODULO as usize + GF_2_16::MAX_EXP as usize + 1];

        // phi : the almost inverse of psi
        // **Notice** phi[max_exp] should be valid
        //
        // phi[p: reduced poly] = q
        // such that q % P = p
        // Therefore, phi[0] = 0.
        let mut phi: Vec<u16> = vec![0; GF_2_16::ORDER as usize];

        // p = the polynomial 1 = \alpha^0
        let mut p = None;

        // 毎回 \alpha^i % P を計算するのが単純な方法だが
        // 次数が大きい場合の除算をするよりも
        // \alpha^{i+1} % P = ((\alpha^i % P) * \alpha) % P
        // の等式を用いて、直前の結果 \alpha^i % P を用いるようにする
        for i in 0u16..=GF_2_16::MAX_EXP {
            if let Some(p_) = p {
                // i > 0
                // p := p * \alpha
                p = Some(p_ * Poly::<GF_2>::from_mono(1, GF_2::ONE));
            } else {
                // i == 0
                // p = the polynomoal 1
                p = Some(Poly::<GF_2>::from_mono(0, GF_2::ONE))
            }

            // v is the binary representetion of the reduced poly of `p`
            let reduced = &(p.unwrap()) % &ppoly;
            p = Some(reduced.clone());

            let rep: GF_2_16_Val = reduced.into();
            let bin_rep: u16 = rep.into();

            // psi(\alpha^i) = \alpha^i % P
            psi[i as usize] = rep;
            psi[i as usize + GF_2_16::MODULO as usize] = rep;

            // phi(\alpha^i % P) = \alpha^i
            phi[bin_rep as usize] = i;
        }

        GF_2_16 { ppoly, psi, phi }
    }

    pub fn ppoly(&self) -> &Poly<GF_2> {
        &self.ppoly
    }

    /// For p, q: reduced poly
    /// p * q is a reduced poly
    pub fn mul(&self, p: GF_2_16_Val, q: GF_2_16_Val) -> GF_2_16_Val {
        let p = u16::from(p);
        let q = u16::from(q);

        // multiplication with the zero equals to the zero
        if p == 0 || q == 0 {
            return 0.into();
        }

        let i: u32 = self.phi[p as usize].into(); // a^i % P = p
        let j: u32 = self.phi[q as usize].into(); // a^j % P = q

        // alpha^i * alpha^j = alpha^{i + j} = \alpha^{(i + j) % modulo}
        // Since \alpha^modulo = 1
        self.psi[(i + j) as usize]
    }

    pub fn mul_inv(&self, p: GF_2_16_Val) -> GF_2_16_Val {
        let p = u16::from(p);

        // pは零元でないと仮定して
        debug_assert!(p != 0);

        // alpha^i = p
        let i = self.phi[p as usize];
        // since a^modulo = 1
        let inv = GF_2_16::MODULO - i;
        self.psi[inv as usize]
    }

    /// For p, q: reduced poly
    /// p + q is obtained by the componentwize addition
    /// Especially, on GF(2), + = XOR
    pub fn add(&self, p: GF_2_16_Val, q: GF_2_16_Val) -> GF_2_16_Val {
        (u16::from(p) ^ u16::from(q)).into()
    }

    /// For p: redureced poly,
    /// -p is the unique value such that p + (-p) = 0
    /// on GF(2), since -v = v, -p=p.
    pub fn add_inv(&self, p: GF_2_16_Val) -> GF_2_16_Val {
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_inv() {
        check_add_inv::<GF_2>();
        check_add_inv::<GF_2_8>();
        check_add_inv::<GF_2_16_Val>();
    }

    #[test]
    fn test_mul_inv() {
        check_mul_inv::<GF_2>();
        check_mul_inv::<GF_2_8>();
        check_mul_inv::<GF_2_16_Val>();
    }

    #[test]
    fn test_distributive_law() {
        check_distributive_law::<GF_2>();
        check_distributive_law::<GF_2_8>();
    }

    fn check_add_inv<F: FiniteField>() {
        for e in F::enumerate() {
            assert_eq!(e + (-e), F::ZERO);
        }
    }

    fn check_mul_inv<F: FiniteField>() {
        for e in F::enumerate() {
            if e != F::ZERO {
                let inv = (&e).mul_inv();
                assert_eq!(e * inv, F::ONE);
            }
        }
    }

    fn check_distributive_law<F: FiniteField>() {
        for e1 in F::enumerate() {
            for e2 in F::enumerate() {
                for e3 in F::enumerate() {
                    assert_eq!(e1 * (e2 + e3), e1 * e2 + e1 * e3);
                }
            }
        }
    }
}
