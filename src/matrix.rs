use crate::field::*;
use std::ops::{Add, Mul};

/*
 * 列ベクトル
 */
pub struct RowV<F: Field> {
    elems: usize,
    inner: Vec<F>,
}

impl<F: Field> RowV<F> {
    pub fn new(elems: usize) -> RowV<F> {
        RowV {
            elems,
            inner: vec![F::ZERO; elems],
        }
    }

    ///
    /// # Safety
    ///
    /// get the reference of the element without index checking
    pub unsafe fn get(&self, idx: usize) -> &F {
        self.inner.get_unchecked(idx)
    }

    ///
    /// # Safety
    ///
    /// get the mutable reference of the element without index checking
    pub unsafe fn get_mut(&mut self, idx: usize) -> &mut F {
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
impl<F: Field> Add for &RowV<F> {
    type Output = RowV<F>;

    fn add(self, rhs: &RowV<F>) -> RowV<F> {
        assert!(self.elems == rhs.elems);

        let elems = self.elems();
        let mut v = vec![F::ZERO; elems];
        for (i, item) in v.iter_mut().enumerate() {
            unsafe { *item = *self.get(i) + *rhs.get(i) }
        }

        RowV { elems, inner: v }
    }
}

/*
 (a1 a2 ... aN) * k = (k*a1 k*a2 ... k*aN)
*/
impl<F: Field> Mul<F> for RowV<F> {
    type Output = RowV<F>;

    fn mul(self, rhs: F) -> RowV<F> {
        let v: Vec<F> = self.iter().map(|e| rhs * *e).collect();
        RowV {
            elems: self.elems(),
            inner: v,
        }
    }
}

/*
 (a1 a2 ... aN) * (b1 b2 ... bN) = (a1 b1) + (a2 b2) + ... + (aN bN)
*/
impl<F: Field> Mul for RowV<F> {
    type Output = F;

    fn mul(self, rhs: RowV<F>) -> F {
        assert!(self.elems() == rhs.elems());
        let mut v = F::ZERO;
        for i in 0..self.elems() {
            unsafe { v = v + (*self.get(i) * *rhs.get(i)) }
        }
        v
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
    cols: usize, // width
    rows: usize, // height

    // 列ベクトルを横に並べているイメージ
    // Vec<RowV<F>>
    inner: Vec<Vec<F>>,
}

impl<F: Field> Matrix<F> {
    /*
     * Make the zero matrix
     */
    pub fn new(cols: usize, rows: usize) -> Matrix<F> {
        let v = vec![vec![F::ZERO; rows]; cols];

        Matrix {
            cols,
            rows,
            inner: v,
        }
    }

    pub fn get(&self, x: usize, y: usize) -> Option<F> {
        if x < self.cols && y < self.rows {
            let row = &self.inner[x];
            Some(row[y])
        } else {
            None
        }
    }

    ///
    /// # Safety
    ///
    /// get the reference of the element without index checking
    pub unsafe fn get_unchecked(&self, x: usize, y: usize) -> &F {
        self.inner.get_unchecked(x).get_unchecked(y)
    }

    ///
    /// # Safety
    ///
    /// get the mutable reference of the element without index checking
    pub unsafe fn get_unchecked_mut(&mut self, x: usize, y: usize) -> &mut F {
        self.inner.get_unchecked_mut(x).get_unchecked_mut(y)
    }

    pub fn get_mut(&mut self, x: usize, y: usize) -> Option<&mut F> {
        if x < self.cols && y < self.rows {
            let row = &mut self.inner[x];
            Some(&mut row[y])
        } else {
            None
        }
    }

    pub fn row_vec(&self, x: usize) -> RowV<F> {
        RowV {
            elems: self.rows,
            inner: self.inner[x].clone(),
        }
    }

    pub fn column_vec(&self, y: usize) -> RowV<F> {
        let mut v = vec![F::ZERO; self.cols];

        for (x, item) in v.iter_mut().enumerate() {
            unsafe {
                *item = *self.get_unchecked(x, y);
            }
        }

        RowV {
            elems: self.cols,
            inner: v,
        }
    }
}

impl<F: Field> Mul for &Matrix<F> {
    type Output = Matrix<F>;

    fn mul(self, rhs: Self) -> Matrix<F> {
        assert!(self.cols == rhs.rows);

        let mut m = Matrix::new(self.cols, rhs.rows);

        for i in 0..self.rows {
            for j in 0..rhs.cols {
                let col_vec = self.column_vec(i);
                let row_vec = rhs.row_vec(j);
                unsafe {
                    *m.get_unchecked_mut(i, j) = col_vec * row_vec;
                }
            }
        }

        m
    }
}
