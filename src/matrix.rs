use crate::field::*;
use crate::vecteur::*;
use std::ops::{Index, IndexMut, Mul};

/*
(h, w)-Matrix is one of the form

( a_00 a_01 a_02             ...     a_0{w-1} )
( a_10 a_11 a_12             ...     a_1{w-1} )
(                            ...              )
( a_{h-1}0 a_{h-1}1 a_{h-1}2 ... a_{h-1}{w-1} )

for 0 <= i <= h-1, 0 <= j <= w-1,
M.at(i, j) equals a_ij
 */
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix<F: Field> {
    // 行ベクトルを縦に並べているイメージ
    inner: Vec<Vecteur<F>>,
}

impl<F: Field> Index<usize> for Matrix<F> {
    type Output = Vecteur<F>;

    fn index(&self, idx: usize) -> &Vecteur<F> {
        &self.inner[idx]
    }
}

impl<F: Field> IndexMut<usize> for Matrix<F> {
    fn index_mut(&mut self, idx: usize) -> &mut Vecteur<F> {
        &mut self.inner[idx]
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct MatrixSize {
    pub height: usize,
    pub width: usize,
}

impl<F: Field> Matrix<F> {
    pub fn iter(&self) -> impl Iterator<Item = &Vecteur<F>> {
        self.inner.iter()
    }

    /*
     * x-行ベクトルとy-行ベクトルを入れ替える
     */
    pub fn swap_column(&mut self, x: usize, y: usize) {
        self.inner.swap(x, y);
    }

    /// ix in ixs ==> 0 <= ix < self.height
    pub fn drop_columns(&mut self, mut ixs: Vec<usize>) {
        ixs.sort();

        // 行列の下の方の行から削除していく
        for ix in ixs.iter().rev() {
            self.inner.remove(*ix);
        }
    }

    /*
     * 行列を入力として
     * その行列が正則行列の場合 ==> 逆行列を Some で包んで返す。
     * 正則行列でない場合 ==> None を返す。
     */
    pub fn inverse(&mut self) -> Option<Matrix<F>> {
        assert!(self.width() == self.height());

        let mut m = Self::identity(self.height());

        for i in 0..self.height() {
            // 軸変換
            if self[i][i] == F::ZERO {
                // この対角成分では掃き出しできないので
                // 他の非零値を下降方向で探す

                let mut found: bool = false;
                for y in i + 1..self.height() {
                    if self[y][i] != F::ZERO {
                        // i-行ベクトルと y-行ベクトルを入れ替える
                        self.swap_column(i, y);
                        m.swap_column(i, y);
                        found = true;
                    }
                }

                // 交換可能な相手がいない場合は
                // この行列が正則行列でないことになるので
                // Noneを返す。
                if !found {
                    return None;
                }
            }

            // 正規化
            let k = self[i][i];
            assert!(k != F::ZERO);
            self[i] = &self[i] * k.mul_inv();
            m[i] = &m[i] * k.mul_inv();

            // 掃き出し
            for x in 0..self.height() {
                if i != x {
                    let k = self[x][i];
                    self[x] = &self[x] - &(&self[i] * k);
                    m[x] = &m[x] - &(&m[i] * k);
                }
            }
        }

        Some(m)
    }

    pub fn height(&self) -> usize {
        self.inner.len()
    }
    pub fn width(&self) -> usize {
        self.inner[0].len()
    }

    /*
     * Make the zero matrix with the given size
     */
    pub fn new(size: MatrixSize) -> Matrix<F> {
        let v = vec![Vecteur::new(size.width); size.height];

        Matrix { inner: v }
    }

    /*
     * Make the identity matrix with the argument
     */
    pub fn identity(size: usize) -> Matrix<F> {
        let mut m = Matrix::new(MatrixSize {
            height: size,
            width: size,
        });

        for i in 0..size {
            m[i][i] = F::ONE;
        }

        m
    }

    /*
     * M.get(x, y) は x番目の行ベクトルのy番目の値
     * 例.
     *   1           2           3           4
     *   5           6           7           8
     *   9          10          11          12
     * なる(3, 4)行列 Aについて
     * A(2, 4)は8である
     */
    pub fn get(&self, i: usize, j: usize) -> Option<F> {
        if i < self.height() && j < self.width() {
            Some(self.inner[i][j])
        } else {
            None
        }
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> Option<&mut F> {
        if i < self.height() && j < self.width() {
            Some(&mut self.inner[i][j])
        } else {
            None
        }
    }

    ///
    /// # Safety
    ///
    /// get the reference of the element without index checking
    pub unsafe fn get_unchecked(&self, i: usize, j: usize) -> &F {
        &self[i][j]
    }

    ///
    /// # Safety
    ///
    /// get the mutable reference of the element without index checking
    pub unsafe fn get_unchecked_mut(&mut self, i: usize, j: usize) -> &mut F {
        &mut self[i][j]
    }

    pub fn column_vec(&self, i: usize) -> &Vecteur<F> {
        &self[i]
    }

    pub fn row_vec(&self, j: usize) -> Vecteur<F> {
        let mut v = vec![F::ZERO; self.height()];

        for (x, item) in v.iter_mut().enumerate() {
            *item = self[x][j];
        }

        Vecteur::from_vec(v)
    }
}

impl<F: Field> Mul for &Matrix<F> {
    type Output = Matrix<F>;

    fn mul(self, rhs: Self) -> Matrix<F> {
        assert!(self.width() == rhs.height());

        let mut m: Matrix<F> = Matrix::new(MatrixSize {
            height: self.height(),
            width: rhs.width(),
        });

        for i in 0..self.height() {
            for j in 0..rhs.width() {
                let col_vec = self.column_vec(i);
                let row_vec = rhs.row_vec(j);
                m[i][j] = col_vec * &row_vec;
            }
        }

        m
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::fin_field::*;

    type MGF2 = Matrix<GF_2>;

    /*
     * 0 1 0 1
     * 1 0 1 0
     * 1 1 1 1
     */
    fn test_matrix1() -> MGF2 {
        let o = GF_2::ZERO;
        let l = GF_2::ONE;
        let mut m = Matrix::new(MatrixSize {
            height: 3,
            width: 4,
        });

        m[0][0] = o;
        m[0][1] = l;
        m[0][2] = o;
        m[0][3] = l;
        m[1][0] = l;
        m[1][1] = o;
        m[1][2] = l;
        m[1][3] = o;
        m[2][0] = l;
        m[2][1] = l;
        m[2][2] = l;
        m[2][3] = l;

        return m;
    }

    /*
     * 0 1 0
     * 1 0 0
     * 0 0 1
     */
    fn test_matrix2() -> MGF2 {
        let o = GF_2::ZERO;
        let l = GF_2::ONE;
        let mut m = Matrix::new(MatrixSize {
            height: 3,
            width: 3,
        });

        m[0][0] = o;
        m[0][1] = l;
        m[0][2] = o;
        m[1][0] = l;
        m[1][1] = o;
        m[1][2] = o;
        m[2][0] = o;
        m[2][1] = o;
        m[2][2] = l;

        return m;
    }

    /*
     * 0 0 1
     * 0 1 0
     * 1 0 0
     */
    fn test_matrix3() -> MGF2 {
        let o = GF_2::ZERO;
        let l = GF_2::ONE;
        let mut m = Matrix::new(MatrixSize {
            height: 3,
            width: 3,
        });

        m[0][0] = o;
        m[0][1] = o;
        m[0][2] = l;
        m[1][0] = o;
        m[1][1] = l;
        m[1][2] = o;
        m[2][0] = l;
        m[2][1] = o;
        m[2][2] = o;

        return m;
    }

    /*
     * 0 1 0
     * 0 0 1
     * 1 0 0
     */
    fn test_matrix4() -> MGF2 {
        let o = GF_2::ZERO;
        let l = GF_2::ONE;
        let mut m = Matrix::new(MatrixSize {
            height: 3,
            width: 3,
        });

        m[0][0] = o;
        m[0][1] = l;
        m[0][2] = o;
        m[1][0] = o;
        m[1][1] = o;
        m[1][2] = l;
        m[2][0] = l;
        m[2][1] = o;
        m[2][2] = o;

        return m;
    }

    /*
     * 0 0 1
     * 1 0 0
     * 0 1 0
     */
    fn test_matrix5() -> MGF2 {
        let o = GF_2::ZERO;
        let l = GF_2::ONE;
        let mut m = Matrix::new(MatrixSize {
            height: 3,
            width: 3,
        });

        m[0][0] = o;
        m[0][1] = o;
        m[0][2] = l;
        m[1][0] = l;
        m[1][1] = o;
        m[1][2] = o;
        m[2][0] = o;
        m[2][1] = l;
        m[2][2] = o;

        return m;
    }

    /*
     * 1 0 0
     * 0 0 1
     * 0 1 0
     */
    fn test_matrix6() -> MGF2 {
        let o = GF_2::ZERO;
        let l = GF_2::ONE;
        let mut m = Matrix::new(MatrixSize {
            height: 3,
            width: 3,
        });

        m[0][0] = l;
        m[0][1] = o;
        m[0][2] = o;
        m[1][0] = o;
        m[1][1] = o;
        m[1][2] = l;
        m[2][0] = o;
        m[2][1] = l;
        m[2][2] = o;

        return m;
    }

    #[test]
    fn height_test() {
        let m = test_matrix1();
        assert_eq!(m.height(), 3);
        assert_eq!(MGF2::identity(4).height(), 4);
    }

    #[test]
    fn width_test() {
        let m = test_matrix1();
        assert_eq!(m.width(), 4);
        assert_eq!(MGF2::identity(4).width(), 4);
    }

    fn make_by_swap<F: Field>(m: &Matrix<F>, i: usize, j: usize) -> Matrix<F> {
        let mut m = m.clone();
        m.swap_column(i, j);
        m
    }

    #[test]
    fn swap_column_test() {
        let i = MGF2::identity(3);
        let a = test_matrix2();
        let b = test_matrix3();

        assert_eq!(make_by_swap(&i, 0, 1), a);
        assert_eq!(make_by_swap(&i, 1, 0), a);
        assert_eq!(make_by_swap(&a, 0, 1), i);
        assert_eq!(make_by_swap(&a, 1, 0), i);

        assert_eq!(make_by_swap(&i, 0, 2), b);
        assert_eq!(make_by_swap(&i, 2, 0), b);
        assert_eq!(make_by_swap(&b, 0, 2), i);
        assert_eq!(make_by_swap(&b, 2, 0), i);
    }

    #[test]
    fn column_vec_test() {
        let m = MGF2::identity(3);
        let v1 = Vecteur::from_vec(vec![GF_2::ONE, GF_2::ZERO, GF_2::ZERO]);
        let v2 = Vecteur::from_vec(vec![GF_2::ZERO, GF_2::ONE, GF_2::ZERO]);
        let v3 = Vecteur::from_vec(vec![GF_2::ZERO, GF_2::ZERO, GF_2::ONE]);

        assert_eq!(m.column_vec(0), &v1);
        assert_eq!(m.column_vec(1), &v2);
        assert_eq!(m.column_vec(2), &v3);
    }

    #[test]
    fn row_vec_test() {
        let m = MGF2::identity(3);
        let v1 = Vecteur::from_vec(vec![GF_2::ONE, GF_2::ZERO, GF_2::ZERO]);
        let v2 = Vecteur::from_vec(vec![GF_2::ZERO, GF_2::ONE, GF_2::ZERO]);
        let v3 = Vecteur::from_vec(vec![GF_2::ZERO, GF_2::ZERO, GF_2::ONE]);

        assert_eq!(m.row_vec(0), v1);
        assert_eq!(m.row_vec(1), v2);
        assert_eq!(m.row_vec(2), v3);
    }

    #[test]
    fn matrix_product_test() {
        let i = MGF2::identity(3);
        let m1 = test_matrix2();
        let m2 = test_matrix3();
        let m3 = test_matrix4();
        let m4 = test_matrix5();
        let m5 = test_matrix6();

        let targets = vec![
            (&i, &i, &i),
            (&m1, &i, &m1),
            (&i, &m1, &m1),
            (&m2, &i, &m2),
            (&i, &m2, &m2),
            (&m3, &i, &m3),
            (&i, &m3, &m3),
            (&m4, &i, &m4),
            (&i, &m4, &m4),
            (&m1, &m2, &m3),
            (&m2, &m1, &m4),
            (&m2, &m3, &m5),
            (&m3, &m2, &m1),
        ];

        for (a, b, c) in targets {
            assert_eq!(a * b, *c);
        }
    }

    #[test]
    fn inverse_matrix_test() {
        let i = MGF2::identity(3);
        let m1 = test_matrix2();
        let m2 = test_matrix3();
        let m3 = test_matrix4();
        let m4 = test_matrix5();
        let m5 = test_matrix6();

        let targets = vec![
            (&i, &i),
            (&m1, &m1),
            (&m2, &m2),
            (&m3, &m4),
            (&m4, &m3),
            (&m5, &m5),
        ];

        for (a, b) in targets {
            assert!(a.clone().inverse().is_some(), "{:?}", a);
            assert_eq!(a.clone().inverse().unwrap(), *b, "\n orig = {:?}", a);
        }
    }
}
