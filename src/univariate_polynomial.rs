use crate::field::*;
use std::collections::BTreeMap;
use std::ops::{Add, Sub, Mul, Div, Rem};

/*
 * 1変数多項式を定義する。
 * 変数はXで表しているとする。
 * 係数は体とする。
 */

#[derive(Clone, Debug)]
pub struct Poly<F: Field> {
    // 多項式の次数をkeyとして係数をvalueとして表現する。
    // X^3 + X + 1 は
    // p[3] = 1, p[1] = 1, p[0] = 1
    // p[3] = 1, p[2] = 0, p[1] = 1, p[0] = 1
    // p[5] = 0, p[3] = 1, p[1] = 1, p[0] = 1
    // など複数の表現を持つ。
    inner: BTreeMap<u32, F>,
}

impl<F: Field> PartialEq for Poly<F> {
    fn eq(&self, rhs: &Self) -> bool {
        self.normalize().inner == rhs.normalize().inner
    }
}

impl<F: Field> Eq for Poly<F> {}

impl<F: Field> Poly<F> {
    fn add_inv(&self) -> Self {
        let mut p = BTreeMap::new();

        for (deg, coeff) in self.inner.iter() {
            p.insert(*deg, coeff.add_inv());
        }

        Poly { inner: p }
    }
    
    fn at(&self, degree: &u32) -> F {
        self.inner.get(degree).copied().unwrap_or(F::ZERO)
    }
    
    pub fn zero() -> Self {
        Poly { inner: BTreeMap::new() }
    }

    pub fn one() -> Self {
        let mut p = BTreeMap::new();
        p.insert(0, F::ONE);
        Poly { inner: p }
    }
    
    fn new() -> Self {
        Poly { inner: BTreeMap::new() }
    }
    
    fn normalize(&self) -> Self {
        let mut p = BTreeMap::new();

        for (d, coeff) in self.inner.iter() {
            if *coeff  != F::ZERO {
                p.insert(*d, *coeff);
            }
        }

        Poly { inner: p }
    }
    
    // 空の多項式 または 係数が全て零 <=> 零元
    fn is_zero(&self) -> bool {
        self.inner.values().all(|coeff| *coeff == F::ZERO)
    }

    // 係数が零でない次数のうち最大のものを返す
    fn degree(&self) -> Option<u32> {
        self.max_mono().map(|(deg, _)| deg)
    }

    // 最大次数を持つ非零項を返す
    fn max_mono(&self) -> Option<(u32, F)> {
        self.inner.iter()
            .rfind(|(_, coeff)| **coeff != F::ZERO)
            .map(|(deg, coeff)| (*deg, *coeff))
    }
    
    pub fn from_btree(b: BTreeMap<u32, F>) -> Self {
        Poly { inner: b }
    }

    pub fn from_vec(v: Vec<(u32, F)>) -> Self {
        let mut b = BTreeMap::new();
        for (k, v) in v {
            b.insert(k, v);
        }
        Poly { inner: b }
    }

    fn from_mono(deg: u32, coeff: F) -> Self {
        let mut b = BTreeMap::new();
        b.insert(deg, coeff);
        Poly { inner: b }
    }
    
    /*
     * long_div(p, q) = (a, b) s.t. p = aq + b where b = 0 or deg(a) < deg(q)
     */
    pub fn long_div(&self, divisor: Self) -> (Poly<F>, Poly<F>) {
        // divisorのdegreeに何か掛けて
        // selfのdegreeを打ち消そうとする

        let mut quotient = Poly::new();
        let mut divided = self.clone();

        loop {
            if divided.is_zero() {
                return (quotient, divided);
            }
            
            if divided.degree() < divisor.degree() {
                return (quotient, divided);
            }

            // divided = c1 x^{d1} + ...
            let (d1, c1) = divided.max_mono().unwrap();

            // divisor = c2 x^{d2} + ...
            let (d2, c2) = divisor.max_mono().unwrap();

            let d = d1 - d2;
            let c = c1 / c2;
            let p = Poly::from_mono(d, c);
            
            // (c2 x^{d2} + ....) * p = c x^{d} = c1 x^{d1} + ...'
            quotient = quotient.clone() + p.clone();
            divided = divided - divisor.clone() * p;
        }
    }
}

impl<F: Field> Add for Poly<F> {
    type Output = Poly<F>;
    
    fn add(self, rhs: Self) -> Poly<F> {
        use std::cmp::max;
        
        let d1 = self.degree();
        let d2 = rhs.degree();
        let mut p = Poly::<F>::new();
        
        if let Some(_d1) = d1 {
            if let Some(_d2) = d2 {
                let d = max(_d1, _d2);
                for i in 0..=d {
                    p.inner.insert(i, self.at(&i) + rhs.at(&i));
                }
                p
            } else { // lhs + 0
                self
            }
        } else { // 0 + rhs
            rhs
        }
    }
}

impl<F: Field> Sub for Poly<F> {
    type Output = Poly<F>;

    // p - q = p + (-q)
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Self) -> Poly<F> {
        self + rhs.add_inv()
    }
}

impl<F: Field> Mul for Poly<F> {
    type Output = Poly<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: Self) -> Poly<F> {
        let mut p: BTreeMap<u32, F> = BTreeMap::new();
        
        for (deg1, coeff1) in self.inner.iter() {
            for (deg2, coeff2) in rhs.inner.iter() {
                let d = deg1 + deg2;
                let c = *coeff1 * *coeff2;

                if let Some(x) = p.get_mut(&d) {
                    *x = *x + c;
                } else {
                    p.insert(d, c);
                }
            }
        }

        Poly { inner: p }
    }
}

impl<F: Field> Div for Poly<F> {
    type Output = Poly<F>;

    fn div(self, rhs: Self) -> Poly<F> {
        let (quotient, _) = self.long_div(rhs);
        quotient
    }
}

impl<F: Field> Rem for Poly<F> {
    type Output = Poly<F>;

    fn rem(self, rhs: Self) -> Poly<F> {
        let (_, reminder) = self.long_div(rhs);
        reminder
    }
}

#[cfg(test)]
mod test {
    use crate::fin_field::*;
    use super::*;
    
    type PGF2 = Poly<GF_2>;

    // 0 x^4 + 0 x^3 + 0 x^2 + 0 x + 0
    fn zero1() -> Poly<GF_2> {
        Poly::<GF_2>::from_vec(
            vec![(4, GF_2::ZERO), (3, GF_2::ZERO), (2, GF_2::ZERO), (1, GF_2::ZERO), (0, GF_2::ZERO)]
        )
    }
    
    // x^3 + x + 1
    fn p1() -> Poly<GF_2> {
        Poly::<GF_2>::from_vec(
            vec![(3, GF_2::ONE), (1, GF_2::ONE), (0, GF_2::ONE)]
        )
    }
    // x^3 + x + 1 (via an other form from p1)
    fn p2() -> Poly<GF_2> {
        Poly::<GF_2>::from_vec(
            vec![(0, GF_2::ONE), (3, GF_2::ONE), (1, GF_2::ONE)]
        )
    }

    fn p3() -> Poly<GF_2> {
        Poly::<GF_2>::from_vec(
            vec![(6, GF_2::ONE), (2, GF_2::ONE), (0, GF_2::ONE)]
        )
    }
    
    #[test]
    fn test_degree() {
        assert_eq!(PGF2::zero().degree(), None);
        assert_eq!(zero1().degree(), None);
        assert_eq!(p1().degree().unwrap(), 3);
        assert_eq!(p2().degree().unwrap(), 3);
    }

    #[test]
    fn test_iszero() {
        assert_eq!(PGF2::zero().is_zero(), true);
        assert_eq!(zero1().is_zero(), true);
        assert_eq!(p1().is_zero(), false);
        assert_eq!(p2().is_zero(), false);
    }

    #[test]
    fn test_eq() {
        assert_ne!(PGF2::zero(), p1());
        assert_eq!(PGF2::zero(), zero1());
        assert_eq!(p1(), p2());
    }

    #[test]
    fn test_add() {
        assert_eq!(PGF2::zero() + zero1(), zero1());
        assert_eq!(p1() + zero1(), p1());
        assert_eq!(p2() + zero1(), p2());
    }

    #[test]
    fn test_mul() {
        let v = vec![
            (zero1() , PGF2::one(), zero1()),
            (p1(),  PGF2::one(), p1()),
            (p1(), p1(), p3()),
            (p1(), p2(), p3()),
        ];

        for (x, y, z) in v {
            assert_eq!(x * y, z);
        }
    }

    #[test]
    fn test_div() {
        let v = vec![
            (p1(), PGF2::one(), p1()),
            (p1(), p1(), PGF2::one()),
            (p1(), p2(), PGF2::one()),
            (PGF2::one(), p1(), PGF2::zero()),
            (p3(), p1(), p1())
        ];

        for (x, y, z) in v {
            assert_eq!(x / y, z);
        }
    }

    #[test]
    fn test_rem() {
        let v = vec![
            (p1(), PGF2::one(), PGF2::zero()),
            (p1(), p1(), PGF2::zero()),
            (p1(), p2(), PGF2::zero()),
            (PGF2::one(), p1(), PGF2::one()),
            (p3(), p1(), PGF2::zero()),
        ];

        for (x, y, z) in v {
            assert_eq!(x % y, z);
        }
    }
}