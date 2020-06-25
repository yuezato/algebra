use crate::field::*;

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

    inner: Vec<Vec<F>> // 列ベクトルを横に並べているイメージ
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
            inner: v
        }
    }

    pub fn at(&self, x: usize, y: usize) -> Option<F> {
        if x < self.cols && y < self.rows {
            let row = &self.inner[x];
            Some(row[y])
        } else {
            None
        }
    }

    pub fn at_mut(&mut self, x: usize, y: usize) -> Option<&mut F> {
        if x < self.cols && y < self.rows {
            let row = &mut self.inner[x];
            Some(&mut row[y])
        } else {
            None
        }        
    }
}
