use crate::matrix::*;
// use crate::vandermonde::*;
use crate::fin_field::*;

/*
1. Vec<u8>から行数とbit数を引数にして行列化する関数を作る
2. encode: data_matrix -> <(block_num, Vec<GF>)> とする関数を作る
3. decode: <(block_num, Vec<GF>)> -> data_matrix とする関数を作る
 3.1. decode時には行列の行位置iについてiにだけ1があるパターンで
      計算をskipしてメモリコピーだけで済ませるやつにする
 3.2. この実装だとsystematic vandermonde行列について
      data blockの計算はすっとばせる
*/

/*
 * ここでデータベクトルからデータマトリックスを作ると
 * 巨大なメモリallocationが発生するかもしれないので
 * この段階ではsliceだけ返す。
 */
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

        // m[i] == 0 0 ... 1 0 .. 0 かつ m[i][l] == l
        // なら i に l をコピーとした方が良い
        if let Some(c) = check_copyable(m[i].as_vec()) {
            // これは m[i] = 0 0 0 ... 1 ... 0 0 の形をしていて
            // m[i][c] = 1 ほかは全て 0 であったといっている。
            // この形の行をデータ行列にかける時には c列目 をコピーすれば良い
            // コピーするだけで良い
            coded[i].copy_from_slice(datam[c]);
        } else {
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
                    let x: F = m[i][j] * F::from_bytes(&data[k..(k + F::BYTE_SIZE)]);
                    let z: F;

                    if j == 0 {
                        // coded[0]は未初期化状態なので場合分けが必要になる
                        z = x;
                    } else {
                        z = F::from_bytes(&coded[i][k..(k + F::BYTE_SIZE)]) + x;
                    }

                    for l in 0..F::BYTE_SIZE {
                        coded[i][k + l] = z.to_byte(l);
                    }
                    k += F::BYTE_SIZE;
                }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::*;
    use crate::vandermonde::*;
    use crate::vecteur::*;

    #[test]
    fn conversion_is_iso1() {
        let r = GF_2_16_Val::PRIMITIVE_ROOT;

        let v = vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4)];

        let v_ = data_to_finfield_vec(&finfield_vec_to_data(&v));

        assert_eq!(v, v_);
    }

    #[test]
    fn conversion_is_iso2() {
        let v: Vec<u8> = (0u8..=0xffu8).collect();
        let v_: Vec<GF_2_16_Val> = data_to_finfield_vec(&v);
        let v_ = finfield_vec_to_data(&v_);

        assert_eq!(v, v_);
    }

    #[test]
    fn optimized_mul_equals_to_mul1() {
        let r = GF_2_16_Val::PRIMITIVE_ROOT;

        let v1 = vandermonde(
            MatrixSize {
                height: 4,
                width: 4,
            },
            vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4)],
        )
        .unwrap();

        let datav: Vec<u8> = (0..80).collect();
        let datav_: Vec<&[u8]> = vec_to_quasi_matrix(&datav, 4);

        let mut datam: Matrix<GF_2_16_Val> = Matrix::new(MatrixSize {
            height: datav_.len(),
            width: datav_[0].len() / 2,
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
    }

    #[test]
    fn optimized_mul_equals_to_mul2() {
        let r = GF_2_16_Val::PRIMITIVE_ROOT;

        let v1 = systematic_vandermonde(
            MatrixSize {
                height: 6,
                width: 4,
            },
            vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4), r.exp(5), r.exp(6)],
        )
        .unwrap();

        let datav: Vec<u8> = (0..120).collect();

        // vec_to_quasi_matrix には data block数を渡す
        let datav_: Vec<&[u8]> = vec_to_quasi_matrix(&datav, v1.width());

        let mut datam: Matrix<GF_2_16_Val> = Matrix::new(MatrixSize {
            height: datav_.len(),
            width: datav_[0].len() / GF_2_16_Val::BYTE_SIZE,
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
    }

    #[test]
    fn test_check_copyable() {
        let r = GF_2_16_Val::PRIMITIVE_ROOT;
        let o = GF_2_16_Val::ZERO;
        let l = GF_2_16_Val::ONE;

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

    #[test]
    fn enc_erase_dec_test1() {
        let r = GF_2_16_Val::PRIMITIVE_ROOT;

        let mds = systematic_vandermonde(
            MatrixSize {
                height: 5,
                width: 4,
            },
            vec![r.exp(1), r.exp(2), r.exp(3), r.exp(4), r.exp(5)],
        )
        .unwrap();

        let datav: Vec<u8> = (0..160).collect();
        let datav: Vec<&[u8]> = vec_to_quasi_matrix(&datav, mds.width());

        let mut encoded_datav: Vec<Vec<u8>> = memory_optimized_mul(&mds, &datav);

        encoded_datav.remove(0);

        // データ消失した位置に対応する行ベクトルを削除する
        let mut erased_mds = mds.clone();
        erased_mds.drop_columns(vec![0]);
        let decoder_matrix = erased_mds.inverse().unwrap();

        let decoded_datav: Vec<Vec<u8>> = memory_optimized_mul2(&decoder_matrix, &encoded_datav);

        for i in 0..datav.len() {
            let a: &[u8] = datav[i];
            let b: &Vec<u8> = &decoded_datav[i];
            assert_eq!(a, &b[..]);
        }
    }
}
