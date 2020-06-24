use crate::field::*;
use std::ops::{Add, Sub, Mul, Div};
use std::fmt::Debug;

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
    const ONE:  GF_2 = GF_2 { value: true };
    
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

pub trait FiniteField: Field
{
    fn enumerate() -> Vec<Self>;
}

impl FiniteField for GF_2 {
    fn enumerate() -> Vec<Self> {
        vec![ GF_2::ZERO, GF_2::ONE ]
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

    fn check_add_inv<F: FiniteField>() where
    {
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
