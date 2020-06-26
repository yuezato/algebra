use crate::field::*;
use crate::univariate_polynomial::*;
use std::convert::{From, Into};
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};

pub trait FiniteField: Field {
    fn enumerate() -> Vec<Self>;
}

/*
 * The section of GF(2)
 */

#[derive(Clone, Copy, Debug)]
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

impl PartialEq for GF_2 {
    fn eq(&self, rhs: &Self) -> bool {
        self.value == rhs.value
    }
}
impl Eq for GF_2 {}

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
}

/*
 * The section of GF(2^16)
 */

lazy_static! {
    // use x^16 + x^12 + x^3 + x^1 + 1.
    static ref GF_2_16_IMPL: GF_2_16 = GF_2_16::new(
        Poly::from_vec( vec![
            (16, GF_2::ONE), (12, GF_2::ONE), (3, GF_2::ONE), (1, GF_2::ONE), (0, GF_2::ONE)
        ])
    );
}

#[allow(non_camel_case_types)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct GF_2_16_Val(u16);

impl Field for GF_2_16_Val {
    const ZERO: GF_2_16_Val = GF_2_16_Val(0);
    const ONE: GF_2_16_Val = GF_2_16_Val(1);

    fn mul_inv(&self) -> GF_2_16_Val {
        GF_2_16_IMPL.mul_inv(*self)
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
    const MAX_EXP: u16 = 0xffff - 2;

    /// Since modulo = max_exp + 1,
    /// for any p % `modulo` = 0, \alpha^p = 1.
    const MODULO: u16 = 0xffff - 1;

    // The order of GF(2^16).
    const ORDER: u16 = 0xffff;

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
        let mut psi: Vec<GF_2_16_Val> = vec![0.into(); (GF_2_16::MAX_EXP + 1) as usize];

        // phi : the almost inverse of psi
        // **Notice** phi[max_exp] should be valid
        //
        // phi[p: reduced poly] = q
        // such that q % P = p
        // Therefore, phi[0] = 0.
        let mut phi: Vec<u16> = vec![0; GF_2_16::ORDER as usize];

        for i in 0u16..=GF_2_16::MAX_EXP {
            // p = \alpha^i
            let p = Poly::<GF_2>::from_mono(i as u32, GF_2::ONE);

            // v is the binary representetion of the reduced poly of `p`
            let reduced = &p % &ppoly;
            let rep: GF_2_16_Val = reduced.into();
            let bin_rep: u16 = rep.into();

            // psi(\alpha^i) = \alpha^i % P
            psi[i as usize] = rep;

            // phi(\alpha^i % P) = \alpha^i
            phi[bin_rep as usize] = i;
        }

        GF_2_16 { ppoly, psi, phi }
    }

    pub fn ppoly(self) -> Poly<GF_2> {
        self.ppoly
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

        let i = self.phi[p as usize]; // a^i % P = p
        let j = self.phi[q as usize]; // a^j % P = q

        // alpha^p * alpha^j = alpha^{i + j} = \alpha^{(i + j) % modulo
        // Since \alpha^modulo = 1
        let r = (i + j) % GF_2_16::MODULO;
        self.psi[r as usize]
    }

    pub fn mul_inv(&self, p: GF_2_16_Val) -> GF_2_16_Val {
        let p = u16::from(p);

        // pは零元でないと仮定して
        assert!(p != 0);

        // alpha^i = p
        let i = self.phi[p as usize];
        // since a^modulo = 1
        let inv = (GF_2_16::MODULO - i) % GF_2_16::MODULO;
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
    fn add_inv() {
        check_add_inv::<GF_2>();
    }

    #[test]
    fn mul_inv() {
        check_mul_inv::<GF_2>();
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
}
