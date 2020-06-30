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

/*
 * 行列積を行う関数
 *
 * 1. 行列積を格納するためのメモリ領域以外を確保しない
 * 2. 行列演算ではなくメモリコピーで済ませられる場合にはそのようにする
 */
pub fn memory_optimized_mul<F: FiniteField>(m: &Matrix<F>, datam: Vec<&[u8]>) -> Vec<Vec<u8>> {
    /*
     * m[i] = (0 0 .. 1 .. 0) かつ m[i][i] = 1の時には
     * d[i] の結果をそのまま使えば良い
     */

    let height = datam.len();
    let width = datam[0].len();

    assert!(height >= width);
    assert!(width % F::BYTE_SIZE == 0);

    let id = Matrix::<F>::identity(width);

    // coded[0..height][0..width]
    // coded: Vec< Vec< uN > > にしたいが
    // それだとendian的に困ったことになる
    let mut coded: Vec<Vec<u8>> = Vec::new();

    for i in 0..height {
        coded.push(Vec::with_capacity(width));

        // FIX!: 必要以上に強い
        // m[i] == 0 0 ... 1 0 .. 0 かつ m[i][l] == l
        // なら i に l をコピーとした方が良い
        if m[i] == id[i] {
            // コピーするだけで良い
            coded[i].copy_from_slice(datam[i]);
        } else {
            // 行ベクトル m[i] と データ行列 d の乗算
            // j: データ行列dを縦に走る変数
            for (j, data) in datam.iter().enumerate() {
                // k: データ行列dを横に走る変数
                let mut k = 0;
                loop {
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
                    if k >= width {
                        break;
                    }
                }
            }
        }
    }

    coded
}

#[cfg(test)]
mod tests {
    #[test]
    fn optimized_mul_equals_to_mul() {}
}
