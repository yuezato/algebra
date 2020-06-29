use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};

/*
 * impl<T : Field> ops::Add<T> for T {
 *    type Output = T;
 *
 *   fn add(self, _rhs: T) {
 *       &self.add(&_rhs)
 *   }
 * }
 * こういう定義をすると
 * E0210を理由にエラーになってしまうので
 * traitをextendする形での実装を行う。
 *
 * 本来であれば 参照同士の演算
 * &Self + &Self や &Self * &Self
 * などを要請したいが、型検査器でエラーにしてしまうようなので
 * 値は積極的にCloneする方向を考える。
 * 参照: https://stackoverflow.com/questions/50660911/
 */
pub trait Field:
    Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + Div<Output = Self>
    + Copy
    + Clone
    + Debug
    + Eq
    + PartialEq
    + Sized
{
    /*
     * the element 0 s.t.
     * 0 + x = x + 0 = 0,
     * 0 * x = x * 0 = 0.
     */
    const ZERO: Self;

    /*
     * the element 1 s.t.
     * 1 * x = x * 1 = x.
     */
    const ONE: Self;

    fn mul_inv(&self) -> Self;

    fn exp(&self, exponent: u32) -> Self {
        if exponent == 0 {
            return Self::ONE;
        }

        let v = self.exp(exponent / 2);
        if exponent % 2 == 0 {
            // x^{2n} = x^n * x^n
            v * v
        } else {
            // x^{2n+1} = x^n * x^n * x
            *self * v * v
        }
    }
}
