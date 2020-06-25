use crate::field::*;
use crate::univariate_polynomial::*;
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Sub};

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

    fn add_inv(&self) -> GF_2 {
        *self
    }
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

impl Sub for GF_2 {
    type Output = GF_2;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: GF_2) -> GF_2 {
        self + rhs.add_inv()
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

    fn add_inv(&self) -> Self {
        self.add_inv()
    }
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

fn conv16(p: Poly<GF_2>) -> u16 {
    let mut v = 0;

    for (deg, coef) in p.iter() {
        if *coef == GF_2::ONE {
            v |= 1 << (*deg)
        }
    }

    v
}

#[allow(non_camel_case_types)]
#[allow(non_snake_case)]
pub struct GF_2_16 {
    // primitive polynomial
    ppoly: Poly<GF_2>,

    // psi : 0 \cup \alpha^i \to GF_2[X]^{< 16}
    // psi is an array of the form u16[2^16]
    psi: Vec<u16>,

    // phi = psi^{-1},
    // phii is an array of the form u16[2^16]
    phi: Vec<u16>,
}

impl GF_2_16 {
    pub fn new(ppoly: Poly<GF_2>) -> GF_2_16 {
        let mut psi = vec![0; 2u32.pow(16) as usize];
        let mut phi = vec![0; 2u32.pow(16) as usize];

        for i in 0..(2u32.pow(16) - 1) {
            let p = Poly::<GF_2>::from_mono(i, GF_2::ONE);
            let v = conv16(p % ppoly.clone());
            psi[i as usize] = v;
            phi[v as usize] = i as u16;
        }

        GF_2_16 { ppoly, psi, phi }
    }

    pub fn ppoly(self) -> Poly<GF_2> {
        self.ppoly
    }

    pub fn mul(self, p: u16, q: u16) -> u16 {
        if p == 0 || q == 0 {
            return 0;
        }

        let p: u32 = self.phi[p as usize].into(); // a^pを意味する。零元ではない
        let q: u32 = self.phi[q as usize].into(); // a^qを意味する。零元ではない
        let r = (p + q) % (2u32.pow(16) - 1);
        self.psi[r as usize]
    }

    pub fn mul_inv(self, p: u16) -> u16 {
        // pは零元でないと仮定して

        let p: u32 = self.phi[p as usize].into();
        let inv = (2u32.pow(16) - 1) - p; // since a^{2^16-1} = 1
        self.psi[inv as usize]
    }

    pub fn add(p: u16, q: u16) -> u16 {
        p ^ q
    }

    pub fn add_inv(p: u16) -> u16 {
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
            let inv = (&e).add_inv();
            assert_eq!(e + inv, F::ZERO);
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
