use crate::field::*;
use std::ops::{Add, Index, IndexMut, Mul, Sub};

/*
 * ベクトル
 */
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Vecteur<F: Field> {
    elems: usize,
    inner: Vec<F>,
}

impl<F: Field> Index<usize> for Vecteur<F> {
    type Output = F;

    fn index(&self, idx: usize) -> &F {
        &self.inner[idx]
    }
}

impl<F: Field> IndexMut<usize> for Vecteur<F> {
    fn index_mut(&mut self, idx: usize) -> &mut F {
        &mut self.inner[idx]
    }
}

impl<F: Field> Vecteur<F> {
    pub fn new(elems: usize) -> Vecteur<F> {
        Self {
            elems,
            inner: vec![F::ZERO; elems],
        }
    }

    ///
    /// # Safety
    ///
    /// get the reference of the element without index checking
    pub unsafe fn get_unchecked(&self, idx: usize) -> &F {
        self.inner.get_unchecked(idx)
    }

    ///
    /// # Safety
    ///
    /// get the mutable reference of the element without index checking
    pub unsafe fn get_unchecked_mut(&mut self, idx: usize) -> &mut F {
        self.inner.get_unchecked_mut(idx)
    }

    pub fn elems(&self) -> usize {
        self.elems
    }

    pub fn iter(&self) -> impl Iterator<Item = &F> {
        self.inner.iter()
    }
}

/*
 (a1 a2 ... aN) + (b1 b2 ... bN) = (a1+b1 a2+b2 ... aN+bN)
*/
impl<F: Field> Add for &Vecteur<F> {
    type Output = Vecteur<F>;

    fn add(self, rhs: &Vecteur<F>) -> Vecteur<F> {
        assert!(self.elems == rhs.elems);

        let elems = self.elems();
        let mut v = vec![F::ZERO; elems];
        for (i, item) in v.iter_mut().enumerate() {
            *item = self[i] + rhs[i];
        }

        Vecteur { elems, inner: v }
    }
}
impl<F: Field> Add for Vecteur<F> {
    type Output = Vecteur<F>;

    fn add(self, rhs: Vecteur<F>) -> Vecteur<F> {
        &self + &rhs
    }
}

/*
 v - w = v + (-1 w)
*/
impl<F: Field> Sub for &Vecteur<F> {
    type Output = Vecteur<F>;

    fn sub(self, rhs: &Vecteur<F>) -> Vecteur<F> {
        self + &(rhs * -F::ONE)
    }
}
impl<F: Field> Sub for Vecteur<F> {
    type Output = Vecteur<F>;

    fn sub(self, rhs: Vecteur<F>) -> Vecteur<F> {
        self + rhs * (-F::ONE)
    }
}

/*
(a1 a2 ... aN) * k = (k*a1 k*a2 ... k*aN)
*/
impl<F: Field> Mul<F> for &Vecteur<F> {
    type Output = Vecteur<F>;

    fn mul(self, rhs: F) -> Vecteur<F> {
        let v: Vec<F> = self.iter().map(|e| rhs * *e).collect();
        Vecteur {
            elems: self.elems(),
            inner: v,
        }
    }
}
impl<F: Field> Mul<F> for Vecteur<F> {
    type Output = Vecteur<F>;

    fn mul(self, rhs: F) -> Vecteur<F> {
        &self * rhs
    }
}

/*
 * (a1 a2 ... aN) * (b1 b2 ... bN) = (a1 b1) + (a2 b2) + ... + (aN bN)
 */
impl<F: Field> Mul for &Vecteur<F> {
    type Output = F;

    fn mul(self, rhs: &Vecteur<F>) -> F {
        assert!(self.elems() == rhs.elems());
        let mut v = F::ZERO;
        for i in 0..self.elems() {
            v = v + (self[i] * rhs[i]);
        }
        v
    }
}
impl<F: Field> Mul for Vecteur<F> {
    type Output = F;

    fn mul(self, rhs: Vecteur<F>) -> F {
        &self * &rhs
    }
}

/*
(h, w)-Matrix is one of the form

( a_00 a_01 a_02             ...     a_0{w-1} )
( a_10 a_11 a_12             ...     a_1{w-1} )
(                            ...              )
( a_{h-1}0 a_{h-1}1 a_{h-1}2 ... a_{h-1}{w-1} )

for 0 <= i <= h-1, 0 <= j <= w-1,
M.at(i, j) equals a_ij
*/
pub struct Matrix<F: Field> {
    width: usize,
    height: usize,

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
    height: usize,
    width: usize,
}

impl<F: Field> Matrix<F> {
    /*
     * x-行ベクトルとy-行ベクトルを入れ替える
     */
    pub fn swap(&mut self, x: usize, y: usize) {
        self.inner.swap(x, y);
    }

    /*
     * 行列を入力として
     * その行列が正則行列の場合 ==> 逆行列を Some で包んで返す。
     * 正則行列でない場合 ==> None を返す。
     */
    pub fn inverse(&mut self) -> Option<Matrix<F>> {
        assert!(self.width() == self.height());

        let mut m = Matrix::new(MatrixSize {
            height: self.height(),
            width: self.width(),
        });

        for i in 0..self.height() {
            // 軸変換
            if self[i][i] == F::ZERO {
                // この対角成分では掃き出しできないので
                // 他の非零値を下降方向で探す

                for y in i + 1..self.height() {
                    if self[y][i] != F::ZERO {
                        // i-行ベクトルと y-行ベクトルを入れ替える
                        self.swap(i, y);
                        m.swap(i, y);
                        break;
                    }
                }

                // 交換可能な相手がいない場合は
                // この行列が正則行列でないことになるので
                // Noneを返す。
                return None;
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

    pub fn width(&self) -> usize {
        self.width
    }
    pub fn height(&self) -> usize {
        self.height
    }

    /*
     * Make the zero matrix
     */
    pub fn new(size: MatrixSize) -> Matrix<F> {
        let v = vec![Vecteur::new(size.width); size.height];

        Matrix {
            width: size.width,
            height: size.height,
            inner: v,
        }
    }

    /*
     * M.get(x, y) は x番目の行ベクトルのy番目の値
     * 例.
     *  16           2           3          13
     *   5          11          10           8
     *   9           7           6          12
     *   4          14          15           1
     * なる(4, 4)行列 Aについて
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

        Vecteur {
            elems: self.height(),
            inner: v,
        }
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
