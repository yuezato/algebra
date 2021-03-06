use crate::fin_field::*;
use crate::matrix::*;
use crate::vandermonde::*;
use std::ops::{Index, IndexMut};

pub fn vec_to_quasi_matrix(v: &[u8], height: usize) -> Vec<&[u8]> {
    assert!(v.len() % height == 0);

    let width = v.len() / height;

    let mut m: Vec<&[u8]> = Vec::with_capacity(height);

    for i in 0..height {
        m.push(&v[(i * width)..((i + 1) * width)])
    }

    m
}

pub fn data_to_finfield_vec<F: FiniteField>(v: &[u8]) -> Vec<F> {
    debug_assert!(v.len() % F::BYTE_SIZE == 0);

    let mut fv = Vec::new();

    let mut i = 0;
    loop {
        if i >= v.len() {
            break;
        }
        fv.push(F::from_bytes(&v[i..i + F::BYTE_SIZE]));
        i += F::BYTE_SIZE;
    }
    fv
}

pub fn finfield_vec_to_data<F: FiniteField>(fv: &[F]) -> Vec<u8> {
    let mut v = Vec::new();

    for f in fv {
        for j in 0..F::BYTE_SIZE {
            v.push(f.to_byte(j));
        }
    }

    v
}

fn check_copyable<F: FiniteField>(v: &[F]) -> Option<usize> {
    let mut res = None;

    for (i, f) in v.iter().enumerate() {
        if *f == F::ONE {
            if res.is_some() {
                // 0* 1 0* 1 の場合
                return None;
            }
            res = Some(i); // 0* 1 の場合
        } else if *f != F::ZERO {
            // 0, 1 以外の値が現れた場合
            return None;
        }
    }

    res
}

pub struct ImmutableMatrix<T> {
    height: usize,
    width: usize,
    inner: Vec<T>,
}

impl<T: Clone> ImmutableMatrix<T> {
    pub fn new(initializing_value: T, size: MatrixSize) -> Self {
        let height = size.height;
        let width = size.width;
        let size = height * width;

        ImmutableMatrix {
            height,
            width,
            inner: vec![initializing_value; size],
        }
    }

    pub fn into_vec(self) -> Vec<T> {
        self.inner
    }

    pub fn to_vec(&self) -> &Vec<T> {
        &self.inner
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.height
    }

    pub fn to_nested_vec(&self) -> Vec<&[T]> {
        let mut v: Vec<&[T]> = Vec::new();
        let step: usize = self.width();
        for i in 0..self.height {
            v.push(&self.inner[i * step..(i + 1) * step]);
        }

        v
    }
}

impl<T> Index<usize> for ImmutableMatrix<T> {
    type Output = [T];

    fn index(&self, idx: usize) -> &[T] {
        let begin = self.width * idx;
        let end = self.width * (idx + 1);
        &self.inner[begin..end]
    }
}

impl<T> IndexMut<usize> for ImmutableMatrix<T> {
    fn index_mut(&mut self, idx: usize) -> &mut [T] {
        let begin = self.width * idx;
        let end = self.width * (idx + 1);
        &mut self.inner[begin..end]
    }
}

pub fn memory_optimized_mul4<F: FiniteField>(m: &Matrix<F>, datam: &[Vec<u8>]) -> Vec<u8> {
    let mut v: Vec<&[u8]> = Vec::new();
    for i in datam {
        v.push(&i);
    }
    memory_optimized_mul3(m, &v)
}

pub fn memory_optimized_mul3<F: FiniteField>(m: &Matrix<F>, datam: &[&[u8]]) -> Vec<u8> {
    let width = datam[0].len();

    assert!(m.height() >= m.width());
    assert!(width % F::BYTE_SIZE == 0);

    // データ行列の行ベクトルの要素数
    // すなわち、符号化ベクトルの行ベクトルの要素数
    let elems = width / F::BYTE_SIZE;

    let mut coded: Vec<u8> = vec![0; m.height() * width];

    for i in 0..m.height() {
        if let Some(c) = check_copyable(m[i].as_vec()) {
            coded[i * width..(i + 1) * width].copy_from_slice(datam[c]);
        } else {
            for (j, data) in datam.iter().enumerate() {
                for k in 0..elems {
                    let x: F = m[i][j] * F::get(data, k);
                    let z: F = if j == 0 {
                        x
                    } else {
                        F::get(&coded, i * elems + k) + x
                    };

                    z.put(&mut coded, i * elems + k);
                }
            }
        }
    }

    coded
}

pub fn memory_optimized_mul2<F: FiniteField>(m: &Matrix<F>, datam: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let mut v: Vec<&[u8]> = Vec::new();
    for i in datam {
        v.push(&i);
    }
    memory_optimized_mul(m, &v)
}

/*
 * 行列積を行う関数
 *
 * 1. 行列積を格納するためのメモリ領域以外を確保しない
 * 2. 行列演算ではなくメモリコピーで済ませられる場合にはそのようにする
 */
pub fn memory_optimized_mul<F: FiniteField>(m: &Matrix<F>, datam: &[&[u8]]) -> Vec<Vec<u8>> {
    /*
     * m[i] = (0 0 .. 1 .. 0) かつ m[i][i] = 1の時には
     * d[i] の結果をそのまま使えば良い
     */

    let width = datam[0].len();

    assert!(m.height() >= m.width());
    assert!(width % F::BYTE_SIZE == 0);

    // coded[0..height][0..width]
    // coded: Vec< Vec< uN > > にしたいが
    // それだとendian的に困ったことになる
    let mut coded: Vec<Vec<u8>> = Vec::new();

    for i in 0..m.height() {
        // coded[i] を reserve
        coded.push(vec![0; width]);

        for (j, data) in datam.iter().enumerate() {
            if m[i][j] == F::ZERO {
                // 結果行列に対して何もしなくて良い
            } else if m[i][j] == F::ONE {
                // データ行ベクトルとのXORを取れば良い
                xor_vecs(&mut coded[i], data);
            } else {
                // 最適化ができない場合
                // k: データ行列の行ベクトルdata上を走る変数
                /*
                for k in 0..(width / F::BYTE_SIZE) {
                    mul_and_xor_vecs(m[i][j], k, &mut coded[i], data);
                }
                 */
                F::mul_then_add(m[i][j], &mut coded[i], data);
            }
        }
    }

    coded
}

pub fn mom3<F: FiniteField>(m: &Matrix<F>, table: &MulTable<F>, datam: &[&[u8]]) -> Vec<Vec<u8>> {
    let width = datam[0].len();

    assert!(m.height() >= m.width());
    assert!(width % F::BYTE_SIZE == 0);

    let mut coded: Vec<Vec<u8>> = Vec::new();

    for i in 0..m.height() {
        coded.push(vec![0; width]);

        for (j, data) in datam.iter().enumerate() {
            if m[i][j] == F::ZERO {
            } else if m[i][j] == F::ONE {
                xor_vecs(&mut coded[i], data);
            } else {
                F::mul_then_add2(&table[(i, j)], &mut coded[i], data);
            }
        }
    }

    coded
}

// v1 ^= v2;
pub fn xor_vecs(v1: &mut [u8], v2: &[u8]) {
    debug_assert!(v1.len() == v2.len());

    unsafe {
        let (prefix1, shorts1, suffix1) = v1.align_to_mut::<u128>();
        let (prefix2, shorts2, suffix2) = v2.align_to::<u128>();

        // dbg!(prefix1.len(), shorts1.len(), suffix1.len());
        // dbg!(prefix2.len(), shorts2.len(), suffix2.len());

        if !prefix1.is_empty() || !suffix1.is_empty() || !prefix2.is_empty() || !suffix2.is_empty()
        {
            // slow implementation
            for i in 0..v1.len() {
                v1[i] ^= v2[i];
            }
        } else {
            for i in 0..shorts1.len() {
                shorts1[i] ^= shorts2[i];
            }
        }
    }
}

// v1 ^= k * v2;
pub fn mul_and_xor_vecs<F: FiniteField>(k: F, idx: usize, v1: &mut [u8], v2: &[u8]) {
    let x: F = k * F::get(v2, idx);
    let z: F = F::get(v1, idx) + x;
    z.put(v1, idx);
}

pub fn dot_prod_row_and_matrix<F: FiniteField>(v: &[F], m: &[&[u8]]) -> Vec<u8> {
    debug_assert!(v.len() == m.len());

    // FIX: 0が零元だと仮定してしまっている
    let mut res: Vec<u8> = vec![0; m[0].len()];

    // i = 0の場合は初期化も兼ねてコピーだけ行う

    // i > 0の場合は足し込んでいく
    for i in 0..m.len() {
        F::mul_then_add(v[i], &mut res[..], m[i]);
    }

    res
}

// t (row) <- v (row) * s (matrix)
pub fn dot_prod_row_and_matrix_into<F: FiniteField>(t: &mut [u8], v: &[F], s: &[&[u8]]) {
    debug_assert!(v.len() == s.len());
    debug_assert!(t.len() == s[0].len());

    for i in 0..s.len() {
        if v[i] == F::ZERO {
            // 何もしなくて良い
        } else if v[i] == F::ONE {
            // 係数が1のxorで済む特別な場合
            xor_vecs(t, s[i]);
        } else {
            F::mul_then_add(v[i], t, s[i]);
        }
    }
}

// t (row) <- v (row) * s (matrix)
pub fn dot_prod_row_and_matrix_by_table_into<F: FiniteField>(
    t: &mut [u8],
    v: &[F],
    s: &[&[u8]],
    table: &MulTable<F>,
    idx: usize,
) {
    debug_assert!(v.len() == s.len());
    debug_assert!(t.len() == s[0].len());

    for i in 0..s.len() {
        if v[i] == F::ZERO {
            // 何もしなくて良い
        } else if v[i] == F::ONE {
            // 係数が1のxorで済む特別な場合
            xor_vecs(t, s[i]);
        } else {
            F::mul_then_add2(&table[(idx, i)], t, s[i]);
        }
    }
}

pub fn mom<F: FiniteField>(m: &Matrix<F>, datam: &[&[u8]]) -> ImmutableMatrix<u8> {
    let width = datam[0].len();

    assert!(m.height() >= m.width());
    assert!(width % F::BYTE_SIZE == 0);

    // We assume
    //   F::ZERO == 0^{F::BYTE_SIZE}
    let mut coded: ImmutableMatrix<u8> = ImmutableMatrix::new(
        0u8,
        MatrixSize {
            height: m.height(),
            width,
        },
    );

    for i in 0..m.height() {
        // 行ベクトル m[i] と データ行列 d の乗算
        // j: データ行列を縦に走る変数
        // data: データ行ベクトル
        //
        // m[i][j] *  データ"行"ベクトルの結果を
        // 結果行列[i][k]に足していく
        //
        // データ"列"ベクトルを使わないので
        // 結果行列[i][k]が求まるわけではないことに注意
        for (j, data) in datam.iter().enumerate() {
            if m[i][j] == F::ZERO {
                // 結果行列に対して何もしなくて良い
            } else if m[i][j] == F::ONE {
                // データ行ベクトルとのXORを取れば良い
                xor_vecs(&mut coded[i], data);
            } else {
                // 最適化ができない場合
                // k: データ行列の行ベクトルdata上を走る変数

                /*
                for k in 0..(width / F::BYTE_SIZE) {
                    mul_and_xor_vecs(m[i][j], k, &mut coded[i], data);
                }
                 */
                F::mul_then_add(m[i][j], &mut coded[i], data);
            }
        }
    }

    coded
}

pub fn mom_another<F: FiniteField>(m: &Matrix<F>, datam: &[&[u8]]) -> ImmutableMatrix<u8> {
    let width = datam[0].len();

    assert!(m.height() >= m.width());
    assert!(width % F::BYTE_SIZE == 0);

    // We assume
    //   F::ZERO == 0^{F::BYTE_SIZE}
    let mut coded: ImmutableMatrix<u8> = ImmutableMatrix::new(
        0u8,
        MatrixSize {
            height: m.height(),
            width,
        },
    );

    for i in 0..m.height() {
        dot_prod_row_and_matrix_into(&mut coded[i], m[i].as_vec(), datam);
    }

    coded
}

pub struct MulTable<F: FiniteField> {
    inner: Vec<Vec<F>>,
    size: MatrixSize,
}

impl<T: FiniteField> MulTable<T> {
    pub fn build(m: &Matrix<T>) -> Self {
        let mut table = MulTable::new(m.size());
        table.setup(m);
        table
    }

    pub fn new(size: MatrixSize) -> Self {
        let elems = size.height * size.width;
        let mut inner = Vec::new();

        for _ in 0..elems {
            let mut v = Vec::with_capacity(T::CARDINALITY);
            unsafe {
                v.set_len(T::CARDINALITY);
            }
            inner.push(v);
        }
        Self { inner, size }
    }

    pub fn setup(&mut self, m: &Matrix<T>) {
        for i in 0..m.height() {
            for j in 0..m.width() {
                for f in T::enumerate() {
                    let r: T = m[i][j] * f;
                    self[(i, j)][f.to_usize()] = r;
                }
            }
        }
    }

    pub fn height(&self) -> usize {
        self.size.height
    }

    pub fn width(&self) -> usize {
        self.size.width
    }
}

impl<T: FiniteField> Index<(usize, usize)> for MulTable<T> {
    type Output = Vec<T>;

    fn index(&self, idx: (usize, usize)) -> &Vec<T> {
        &self.inner[idx.0 * self.width() + idx.1]
    }
}

impl<T: FiniteField> IndexMut<(usize, usize)> for MulTable<T> {
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut Vec<T> {
        let pos = idx.0 * self.width() + idx.1;
        &mut self.inner[pos]
    }
}

pub fn mom2<F: FiniteField>(
    m: &Matrix<F>,
    table: &MulTable<F>,
    datam: &[&[u8]],
) -> ImmutableMatrix<u8> {
    let width = datam[0].len();

    debug_assert!(m.height() >= m.width());
    debug_assert!(width % F::BYTE_SIZE == 0);

    let mut coded: ImmutableMatrix<u8> = ImmutableMatrix::new(
        0u8,
        MatrixSize {
            height: m.height(),
            width,
        },
    );

    for i in 0..m.height() {
        for (j, data) in datam.iter().enumerate() {
            if m[i][j] == F::ZERO {
            } else if m[i][j] == F::ONE {
                xor_vecs(&mut coded[i], data);
            } else {
                F::mul_then_add2(&table[(i, j)], &mut coded[i], data);
            }
        }
    }

    coded
}

pub fn matrix_mul<F: FiniteField>(m: &Matrix<F>, datam: &Matrix<F>) -> Vec<Vec<u8>> {
    let width = datam.width();

    assert!(m.height() >= m.width());

    // coded[0..height][0..width]
    // coded: Vec< Vec< uN > > にしたいが
    // それだとendian的に困ったことになる
    let mut coded: Vec<Vec<u8>> = Vec::new();

    for i in 0..m.height() {
        // coded[i] を reserve
        coded.push(vec![0; width * F::BYTE_SIZE]);

        // 行ベクトル m[i] と データ行列 d の乗算
        // j: データ行列dを縦に走る変数
        // data: データ行ベクトル
        for (j, data) in datam.iter().enumerate() {
            // k: データ行列行ベクトルdataを横に走る変数
            let mut k = 0;
            loop {
                if k >= width {
                    break;
                }
                let x: F = m[i][j] * data[k];
                let z: F;

                if j == 0 {
                    // coded[0]は未初期化状態なので場合分けが必要になる
                    z = x;
                } else {
                    let old = F::from_bytes(&coded[i][k * F::BYTE_SIZE..(k + 1) * F::BYTE_SIZE]);
                    z = old + x;
                }

                for l in 0..F::BYTE_SIZE {
                    coded[i][k * F::BYTE_SIZE + l] = z.to_byte(l);
                }

                k += 1;
            }
        }
    }

    coded
}

pub fn matrix_mul2<F: FiniteField>(m: &Matrix<F>, datam: &Matrix<F>) -> Matrix<F> {
    let width = datam.width();

    assert!(m.height() >= m.width());

    // coded[0..height][0..width]
    // coded: Vec< Vec< uN > > にしたいが
    // それだとendian的に困ったことになる
    let mut coded = Matrix::new(MatrixSize {
        height: m.height(),
        width,
    });

    for i in 0..m.height() {
        // 行ベクトル m[i] と データ行列 d の乗算
        // j: データ行列dを縦に走る変数
        // data: データ行ベクトル
        for (j, data) in datam.iter().enumerate() {
            // k: データ行列行ベクトルdataを横に走る変数
            let mut k = 0;
            loop {
                if k >= width {
                    break;
                }
                let x: F = m[i][j] * data[k];
                let z: F;

                if j == 0 {
                    // coded[0]は未初期化状態なので場合分けが必要になる
                    z = x;
                } else {
                    let old = coded[i][k];
                    z = old + x;
                }

                coded[i][k] = z;

                k += 1;
            }
        }
    }

    coded
}

pub struct Encoded(usize, Vec<u8>);

impl Encoded {
    pub fn block_num(&self) -> usize {
        self.0
    }

    pub fn data(&self) -> &[u8] {
        &self.1
    }
}

#[derive(Debug, Clone)]
pub struct Generator<F: FiniteField>(Matrix<F>);

impl<T: FiniteField> Generator<T> {
    pub fn new(m: Matrix<T>) -> Self {
        Generator(m)
    }

    pub fn matrix(&self) -> &Matrix<T> {
        &self.0
    }

    pub fn take_matrix(self) -> Matrix<T> {
        self.0
    }
}

// TODO: 行列のサイズは大したことがないので
// ここで返してしまっても良い気がする。
#[allow(non_snake_case)]
pub fn encode_by_RSV<F: FiniteField + HasPrimitiveElement>(
    data_size: usize,
    parity_size: usize,
    data: &[u8],
) -> (Generator<F>, Vec<Encoded>) {
    // FIX:
    // 本当は data_size + parity_size <= F::CARDINALITY - 2
    // (-2は0と1を取り除くため)を確認する必要あり。

    let velems: Vec<F> = (1..=(data_size + parity_size))
        .map(|i| F::PRIMITIVE_ELEMENT.exp(i as u32))
        .collect();

    let mds = modified_systematic_vandermonde(
        // systematic_vandermonde(
        MatrixSize {
            height: data_size + parity_size,
            width: data_size,
        },
        &velems,
    )
    .unwrap();

    let datav: Vec<&[u8]> = vec_to_quasi_matrix(&data, mds.width());

    let encoded_datav = memory_optimized_mul(&mds, &datav);

    let mut result: Vec<Encoded> = Vec::new();
    for (i, data) in encoded_datav.into_iter().enumerate() {
        result.push(Encoded(i, data));
    }

    (Generator(mds), result)
}

#[allow(non_snake_case)]
pub fn encode_by_RSV_table<F: FiniteField + HasPrimitiveElement>(
    matrix: &Matrix<F>,
    table: &MulTable<F>,
    data: &[u8],
) -> Vec<Encoded> {
    let data_size = matrix.width();
    // FIX:
    // 本当は data_size + parity_size <= F::CARDINALITY - 2
    // (-2は0と1を取り除くため)を確認する必要あり。

    let datav: Vec<&[u8]> = vec_to_quasi_matrix(&data, data_size);

    let encoded_datav = mom3(matrix, table, &datav);

    let mut result: Vec<Encoded> = Vec::new();
    for (i, data) in encoded_datav.into_iter().enumerate() {
        result.push(Encoded(i, data));
    }

    result
}

#[allow(non_snake_case)]
pub fn decode_by_RSV<F: FiniteField>(generator: Generator<F>, data: Vec<Encoded>) -> Vec<u8> {
    let mut generator: Matrix<F> = generator.0;
    let total_block = generator.height();

    let mut erased_columns: Vec<usize> = Vec::new();
    let alive_columns: Vec<usize> = data.iter().map(|e| e.0).collect();
    for i in 0..total_block {
        if !alive_columns.contains(&i) {
            erased_columns.push(i);
        }
    }

    generator.drop_columns(erased_columns);
    let decoder_matrix = generator.inverse().unwrap();
    let encoded_data: Vec<&[u8]> = data.iter().map(|e| &e.1[..]).collect();

    mom(&decoder_matrix, &encoded_data).into_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vecteur::*;

    #[test]
    fn conversion_is_iso1() {
        fn testfunc<F: FiniteField>(r: F) {
            let v = vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4)];
            let v_ = data_to_finfield_vec(&finfield_vec_to_data(&v));
            assert_eq!(v, v_);
        }

        testfunc::<GF_2_8>(GF_2_8::PRIMITIVE_ROOT);
        testfunc::<GF_2_16_Val>(GF_2_16_Val::PRIMITIVE_ROOT);
    }

    #[test]
    fn conversion_is_iso2() {
        fn testfunc<F: FiniteField>() {
            let v: Vec<u8> = (0u8..=0xffu8).collect();
            let v_: Vec<F> = data_to_finfield_vec(&v);
            let v_ = finfield_vec_to_data(&v_);

            assert_eq!(v, v_);
        }

        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();
    }

    #[test]
    fn optimized_mul_equals_to_mul1() {
        testfunc(GF_2_8::PRIMITIVE_ROOT);
        testfunc(GF_2_16_Val::PRIMITIVE_ROOT);

        fn testfunc<F: FiniteField>(r: F) {
            let v1 = vandermonde(
                MatrixSize {
                    height: 4,
                    width: 4,
                },
                &vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4)],
            )
            .unwrap();

            let datav: Vec<u8> = (0..80).collect();
            let datav_: Vec<&[u8]> = vec_to_quasi_matrix(&datav, v1.width());

            let mut datam: Matrix<F> = Matrix::new(MatrixSize {
                height: datav_.len(),
                width: datav_[0].len() / F::BYTE_SIZE,
            });
            for i in 0..datam.height() {
                let v = data_to_finfield_vec(&datav_[i]);
                datam[i] = Vecteur::from_vec(v);
            }

            for i in 0..datam.height() {
                assert_eq!(datam[i].as_vec(), &data_to_finfield_vec(datav_[i]));
            }

            let r1 = matrix_mul2(&v1, &datam);
            let r2 = matrix_mul(&v1, &datam);
            let r3 = memory_optimized_mul(&v1, &datav_);
            let r4 = memory_optimized_mul3(&v1, &datav_);
            let r5 = mom(&v1, &datav_);

            let r6 = {
                let mut table = MulTable::new(v1.size());
                table.setup(&v1);
                mom2(&v1, &table, &datav_)
            };

            let r7 = {
                let mut table = MulTable::new(v1.size());
                table.setup(&v1);
                mom3(&v1, &table, &datav_)
            };

            let r8 = mom_another(&v1, &datav_);

            for i in 0..r1.height() {
                let data = &r2[i];
                assert_eq!(r1[i].as_vec(), &data_to_finfield_vec(data));
            }

            for i in 0..r2.len() {
                assert_eq!(r2[i], r3[i]);
            }

            for i in 0..r1.height() {
                assert_eq!(r1[i].as_vec(), &data_to_finfield_vec(&r3[i]));
                assert_eq!(finfield_vec_to_data(r1[i].as_vec()), r3[i]);
            }

            assert_eq!(r3.concat(), r4);
            assert_eq!(r5.to_vec(), &r4);
            assert_eq!(r6.into_vec(), r4);
            assert_eq!(r7.concat(), r4);
            assert_eq!(r5.to_vec(), r8.to_vec());
        }
    }

    #[test]
    fn optimized_mul_equals_to_mul2() {
        testfunc(GF_2_8::PRIMITIVE_ROOT);
        testfunc(GF_2_16_Val::PRIMITIVE_ROOT);

        fn testfunc<F: FiniteField>(r: F) {
            let v1 = modified_systematic_vandermonde(
                // systematic_vandermonde(
                MatrixSize {
                    height: 6,
                    width: 4,
                },
                &vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4), r.exp(5), r.exp(6)],
            )
            .unwrap();

            let datav: Vec<u8> = (0..120).collect();

            // vec_to_quasi_matrix には data block数を渡す
            let datav_: Vec<&[u8]> = vec_to_quasi_matrix(&datav, v1.width());

            let mut datam: Matrix<F> = Matrix::new(MatrixSize {
                height: datav_.len(),
                width: datav_[0].len() / F::BYTE_SIZE,
            });
            for i in 0..datam.height() {
                let v = data_to_finfield_vec(&datav_[i]);
                datam[i] = Vecteur::from_vec(v);
            }

            for i in 0..datam.height() {
                assert_eq!(datam[i].as_vec(), &data_to_finfield_vec(datav_[i]));
            }

            let r1 = matrix_mul2(&v1, &datam);
            let r2 = matrix_mul(&v1, &datam);
            let r3 = memory_optimized_mul(&v1, &datav_);
            let r4 = memory_optimized_mul3(&v1, &datav_);
            let r5 = mom(&v1, &datav_);
            let r6 = mom_another(&v1, &datav_);

            for i in 0..r1.height() {
                let data = &r2[i];
                assert_eq!(r1[i].as_vec(), &data_to_finfield_vec(data));
            }

            for i in 0..r2.len() {
                assert_eq!(r2[i], r3[i]);
            }

            for i in 0..r1.height() {
                assert_eq!(r1[i].as_vec(), &data_to_finfield_vec(&r3[i]));
                assert_eq!(finfield_vec_to_data(r1[i].as_vec()), r3[i]);
            }

            assert_eq!(r3.concat(), r4);
            assert_eq!(r5.to_vec(), &r4);
            assert_eq!(r5.to_vec(), r6.to_vec());
        }
    }

    #[test]
    fn test_check_copyable() {
        testfunc(GF_2_8::PRIMITIVE_ROOT);
        testfunc(GF_2_16_Val::PRIMITIVE_ROOT);

        fn testfunc<F: FiniteField>(r: F) {
            let o = F::ZERO;
            let l = F::ONE;

            let v1 = vec![o, o];
            assert_eq!(check_copyable(&v1), None);

            let v2 = vec![o, l];
            assert_eq!(check_copyable(&v2), Some(1));

            let v3 = vec![l, o];
            assert_eq!(check_copyable(&v3), Some(0));

            let v4 = vec![l, o, r];
            assert_eq!(check_copyable(&v4), None);

            let v5 = vec![o, o, r];
            assert_eq!(check_copyable(&v5), None);

            let v6 = vec![r, o, l];
            assert_eq!(check_copyable(&v6), None);
        }
    }

    #[test]
    fn enc_erase_dec_test1() {
        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();

        fn testfunc<F: FiniteField + HasPrimitiveElement>() {
            let data_size = 4;
            let parity_size = 1;

            let datav: Vec<u8> = (0..160).collect();
            let (generator, mut encoded) = encode_by_RSV::<F>(data_size, parity_size, &datav);

            encoded.remove(0);

            let decoded: Vec<u8> = decode_by_RSV::<F>(generator, encoded);
            assert_eq!(decoded, datav);
        }
    }
}
