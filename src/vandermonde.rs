use crate::field::*;
use crate::matrix::*;
use std::convert::TryInto;

/*
 * vandermonde(size={height: m, width: n}, v=[a, b, c, ..., x]) |v| = m
 * is
 * (1 a a^2 a^3 a^4 ... a^n)
 * (1 b b^2 b^3 b^4 ... b^n)
 * (1 c c^2 c^3 c^4 ... c^n)
 *           ...
 * (1 x x^n x^3 x^4 ... x^n)
 *
 * Note: n == m is not needed
 */
pub fn vandermonde<F: Field>(size: MatrixSize, v: &[F]) -> Option<Matrix<F>> {
    let mut m = Matrix::new(size);

    // vの要素を縦に並べることになるのでサイズ検査
    if size.height != v.len() {
        return None;
    }

    // vは相異なる要素からなっていなければならない
    for i in 0..v.len() {
        if v[i] == F::ONE {
            return None;
        }
        for j in (i + 1)..v.len() {
            if v[i] == v[j] {
                return None;
            }
        }
    }

    for i in 0..size.height {
        for j in 0..size.width {
            let e: u32 = j.try_into().unwrap();
            m[i][j] = v[i].exp(e);
        }
    }

    Some(m)
}

pub fn systematic_vandermonde<F: Field>(size: MatrixSize, v: &[F]) -> Option<Matrix<F>> {
    let m = vandermonde(size, v);

    if let Some(m) = m {
        let mut sub = m.clone();
        sub.drop_columns((size.width..size.height).collect());
        let inv = sub.inverse().unwrap();
        Some(&m * &inv)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fin_field::*;

    #[test]
    fn test_inverse_vandermonde() {
        let r = GF_2_16_Val::PRIMITIVE_ROOT;

        let v1 = vandermonde(
            MatrixSize {
                height: 4,
                width: 4,
            },
            &vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4)],
        )
        .unwrap();

        let v1_inv = v1.clone().inverse().unwrap();
        assert_eq!(&v1 * &v1_inv, Matrix::identity(4));

        let v2 = vandermonde(
            MatrixSize {
                height: 4,
                width: 4,
            },
            &vec![r.exp(2), r.exp(2), r.exp(3), r.exp(4)],
        );

        assert!(v2.is_none());

        let v3 = vandermonde(
            MatrixSize {
                height: 5,
                width: 4,
            },
            &vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4), r.exp(5)],
        )
        .unwrap();

        // どの1行を取り除いた正方行列も正則であることの確認
        for i in 0..5 {
            let mut v = v3.clone();
            v.drop_columns(vec![i]);

            let v_inv = v.clone().inverse().unwrap();
            assert_eq!(&v * &v_inv, Matrix::identity(4));
        }
    }

    #[test]
    fn systematic_vandermonde_test() {
        let r = GF_2_16_Val::PRIMITIVE_ROOT;

        let mut sv = systematic_vandermonde(
            MatrixSize {
                height: 5,
                width: 4,
            },
            &vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4), r.exp(5)],
        )
        .unwrap();

        // 上の方が単位行列になっていることの確認
        sv.drop_columns(vec![4]);

        assert_eq!(sv, Matrix::identity(4));
    }
}
