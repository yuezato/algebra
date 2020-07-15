use crate::erasure_code::{decode_matrix, AliveBlocks};
use crate::field::*;
use crate::fin_field::*;
use crate::matrix::*;
use crate::reed_solomon::*;
use crate::univariate_polynomial::*;
use fast_array_ops::*;

pub fn xor_vec(v1: &mut [u8], v2: &[u8]) {
    debug_assert!(v1.len() == v2.len());

    /*
    for i in 0..v1.len() {
        v1[i] ^= v2[i];
    }
     */

    unsafe {
        let (prefix1, shorts1, suffix1) = v1.align_to_mut::<u128>();
        let (prefix2, shorts2, suffix2) = v2.align_to::<u128>();

        if !prefix1.is_empty() || !suffix1.is_empty() || !prefix2.is_empty() || !suffix2.is_empty()
        {
            println!("slow implementation");
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

/*
#[target_feature(enable = "avx2")]
pub unsafe fn avx2_xor1(dst: &mut [u8], src: &[u8]) {
    use std::arch::x86_64::*;

    debug_assert!(dst.len() == src.len());
    debug_assert!(dst.len() % 32 == 0);

    let mut di: *mut __m256i = dst.as_mut_ptr() as *mut __m256i;
    let mut si: *const __m256i = src.as_ptr() as *const __m256i;

    for _ in 0..(dst.len() / 32) {
        let d = _mm256_load_si256(di);
        let s = _mm256_load_si256(si);
        let xored = _mm256_xor_si256(d, s);
        _mm256_store_si256(di, xored);
        di = di.add(1);
        si = si.add(1);
    }
}
 */

fn make_sched(v: &[u8]) -> Vec<usize> {
    let mut r = Vec::new();

    for (i, ve) in v.iter().enumerate() {
        for j in 0..8 {
            if (ve >> j) & 1 == 1 {
                r.push(i * 8 + (7 - j));
            }
        }
    }

    r.sort();

    r
}

/*
 * bitmatrixの積はどう実装する??
 *
 * まず行列の方は GF_2_8 上の (d+p, d) 行列をうけとるとして
 * データの方は d*8 になってないといけないのか
 *
 * 元データが [x1, x2, x3, ..., x320] だとすると
 *
 * d=4だとすると縦を32にしたいので
 * x001 x002 x003 ... x010
 * x011 x012 x013 ... x020
 * ..
 * x311 x312 x313 ... x320
 * になるのかな
 *
 * u128 = u8*16 で SIMD したいんだとすると
 * データ行列の横幅は16の倍数にしておきたいのか
 */
pub fn bit_matrix_row_prod(v: &[u8], data: &[&[u8]]) -> Vec<u8> {
    debug_assert!(v.len() * 8 == data.len());

    let mut r: Vec<u8> = vec![0; data[0].len()];
    let sched: Vec<usize> = make_sched(&v);

    for i in sched {
        xor_vec(&mut r, &data[i]);
    }

    r
}

lazy_static! {
    // use x^8 + x^4 + x^3 + x^2 + 1.
    pub static ref BIT_GF_2_8_IMPL: Bit_GF_2_8_impl = Bit_GF_2_8_impl::build(
        Poly::from_vec( vec![
            (8, GF_2::ONE), (4, GF_2::ONE), (3, GF_2::ONE), (2, GF_2::ONE), (0, GF_2::ONE)
        ])
    );
}

pub fn expand(v: &[[u8; 8]]) -> Vec<Vec<u8>> {
    let mut w = Vec::new();

    for j in 0..8 {
        let mut tmp = Vec::new();
        for ve in v {
            tmp.push(ve[j]);
        }
        w.push(tmp);
    }

    w
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Command {
    // Copy(from, to) where
    // 0 <= from < width of data matrix = 8 * width (of a generator),
    // 0 <= to < 8 * height (of a generator).
    Copy(usize, usize),

    // Xor(from, to)
    Xor(usize, usize),
}

impl Command {
    pub fn to(&self) -> usize {
        match self {
            Command::Copy(_, to) => *to,
            Command::Xor(_, to) => *to,
        }
    }
}

pub fn matrix_to_commands(m: &Matrix<GF_2_8>) -> Vec<Command> {
    let mut commands: Vec<Command> = Vec::new();

    for i in 0..m.height() {
        let row: Vec<u8> = m[i].as_vec().iter().map(|e| u8::from(*e)).collect();
        let v: Vec<[u8; 8]> = row.iter().map(|e| BIT_GF_2_8_IMPL.fast_conv(*e)).collect();
        let v: Vec<Vec<u8>> = expand(&v);

        for (j, e) in v.iter().enumerate() {
            commands.append(&mut gen_command(e, i * 8 + j));
        }
    }

    commands
}

pub fn gen_command(v: &[u8], to: usize) -> Vec<Command> {
    let mut r: Vec<Command> = Vec::new();

    let mut to_copy = true;
    for (i, ve) in v.iter().enumerate() {
        for j in 0..8 {
            let from = i * 8 + j;
            if (ve >> (7 - j)) & 1 == 1 {
                if to_copy {
                    to_copy = false;
                    r.push(Command::Copy(from, to));
                } else {
                    r.push(Command::Xor(from, to));
                }
            }
        }
    }

    r
}

pub struct BitEnc(pub usize, Vec<Vec<u8>>);

pub fn bitmatrix_enc1(m: &Matrix<GF_2_8>, data: &[u8]) -> Vec<BitEnc> {
    debug_assert!(data.len() % (m.width() * 8) == 0);

    let mut result: Vec<Vec<u8>> = to_bitmatrix3(m, &vec_to_quasi_matrix(&data, 8 * m.width()));

    let mut output = Vec::new();

    for i in 0..m.height() {
        let mut tmp: Vec<Vec<u8>> = Vec::new();
        for _ in 0..8 {
            tmp.push(result.remove(0));
        }
        output.push(BitEnc(i, tmp));
    }

    output
}

pub fn bitmatrix_enc2(
    commands: Vec<Command>,
    nr_data: usize,
    nr_parity: usize,
    data: &[u8],
) -> Vec<BitEnc> {
    debug_assert!(data.len() % (nr_data * 8) == 0);

    let mut result: Vec<Vec<u8>> = to_bitmatrix4(
        commands,
        nr_data,
        nr_parity,
        &vec_to_quasi_matrix(&data, 8 * nr_data),
    );

    let mut output = Vec::new();

    for i in 0..(nr_data + nr_parity) {
        let mut tmp: Vec<Vec<u8>> = Vec::new();
        for _ in 0..8 {
            tmp.push(result.remove(0));
        }
        output.push(BitEnc(i, tmp));
    }

    output
}

pub fn bitmatrix_enc3(
    commands: Vec<Command>,
    nr_data: usize,
    nr_parity: usize,
    data: &[u8],
) -> Vec<BitEnc> {
    debug_assert!(data.len() % (nr_data * 8) == 0);

    let mut result: Vec<Vec<u8>> = to_bitmatrix4_optim(
        commands,
        nr_data,
        nr_parity,
        &vec_to_quasi_matrix(&data, 8 * nr_data),
    );

    let mut output = Vec::new();

    for i in 0..(nr_data + nr_parity) {
        let mut tmp: Vec<Vec<u8>> = Vec::new();
        for _ in 0..8 {
            tmp.push(result.remove(0));
        }
        output.push(BitEnc(i, tmp));
    }

    output
}

pub fn bitmatrix_dec1(m: &Matrix<GF_2_8>, data: &[BitEnc]) -> Vec<u8> {
    let mut lives: Vec<bool> = vec![false; m.height()];
    for e in data {
        lives[e.0] = true;
    }

    let alive = AliveBlocks::from_boolvec(&lives);

    let invm = decode_matrix(m.clone(), &alive).unwrap();

    let mut encoded_data: Vec<&[u8]> = Vec::new();

    for e in data {
        for f in &e.1 {
            encoded_data.push(f);
        }
    }

    to_bitmatrix3(&invm, &encoded_data).concat()
}

pub fn inv_commands(m: &Matrix<GF_2_8>, data: &[BitEnc]) -> Vec<Command> {
    let mut lives: Vec<bool> = vec![false; m.height()];
    for e in data {
        lives[e.0] = true;
    }

    let alive = AliveBlocks::from_boolvec(&lives);

    let invm = decode_matrix(m.clone(), &alive).unwrap();

    matrix_to_commands(&invm)
}

pub fn bitmatrix_dec2(commands: Vec<Command>, data: &[BitEnc]) -> Vec<u8> {
    let mut encoded_data: Vec<&[u8]> = Vec::new();

    for e in data {
        for f in &e.1 {
            encoded_data.push(f);
        }
    }

    to_bitmatrix5(commands, &encoded_data)
}

pub fn to_bitmatrix5(commands: Vec<Command>, data: &[&[u8]]) -> Vec<u8> {
    let width = data[0].len();
    // (data.len(), width) の行列のつもり
    let mut output: Vec<u8> = Vec::with_capacity(data.len() * width);
    unsafe {
        output.set_len(data.len() * width);
    }

    for c in commands {
        if let Command::Copy(from, to) = c {
            let mslice: &mut [u8] = &mut output[(width * to)..(width * (to + 1))];
            // mslice.copy_from_slice(data[from]);
            fast_array_copy(mslice, data[from]);
        } else if let Command::Xor(from, to) = c {
            let mslice: &mut [u8] = &mut output[(width * to)..(width * (to + 1))];
            // xor_vec(mslice, data[from]);
            fast_array_xor(mslice, data[from]);
        }
    }

    output
}

pub fn bitmatrix_dec3(
    commands: Vec<Command>,
    data: &[BitEnc],
    last_idx_of_erased_data: usize,
) -> Vec<u8> {
    let mut encoded_data: Vec<&[u8]> = Vec::new();
    let topmost_parity_idx = data.len();

    let mut topmost_parities: &[Vec<u8>] = &[];

    for e in data {
        if e.0 == topmost_parity_idx {
            topmost_parities = &e.1;
        }
        for f in &e.1 {
            encoded_data.push(f);
        }
    }

    if topmost_parities.is_empty() {
        panic!("The topmost parity data was lost");
    }

    to_bitmatrix6(
        commands,
        &encoded_data,
        topmost_parities,
        last_idx_of_erased_data,
    )
}

pub fn last_index_of_erased_data(nr_data: usize, survived_data: &[BitEnc]) -> Option<usize> {
    let survived_blocks: Vec<usize> = survived_data.iter().map(|e| e.0).collect();
    (0..nr_data).filter(|i| !survived_blocks.contains(i)).max()
}

// topmost parity が生き残っている場合
// 今の実装だと不味い！！
// データ行列の最後をやるのではなくて
// データ行列の消失してる最後をやる必要がある
pub fn to_bitmatrix6(
    filtered_commands: Vec<Command>,
    data: &[&[u8]],
    parity: &[Vec<u8>],
    last_idx_of_erased_data: usize,
) -> Vec<u8> {
    let width = data[0].len();
    // (data.len(), width) の行列のつもり
    let mut output: Vec<u8> = Vec::with_capacity(data.len() * width);
    unsafe {
        output.set_len(data.len() * width);
    }

    // このcommandsには
    // 元データの最後の1列相当は含まれていないものとする
    for c in filtered_commands {
        if let Command::Copy(from, to) = c {
            let mslice: &mut [u8] = &mut output[(width * to)..(width * (to + 1))];
            // mslice.copy_from_slice(data[from]);
            fast_array_copy(mslice, &data[from]);
        } else if let Command::Xor(from, to) = c {
            let mslice: &mut [u8] = &mut output[(width * to)..(width * (to + 1))];
            fast_array_xor(mslice, data[from]);
        }
    }

    // 最後の一列は data.len()-1 個のデータ行と
    // topmost parity行の足し合わせで復元できる
    let original_width = width * 8;
    let original_height = data.len() / 8;

    let mslice: &mut [u8] = unsafe {
        let ptr: *mut u8 = output.as_ptr() as *mut u8;
        let ptr: *mut u8 = ptr.add(original_width * last_idx_of_erased_data);
        std::slice::from_raw_parts_mut(ptr, original_width)
        // &mut output[original_width*last_idx..original_width*original_height];
    };

    // BitEnc(usize, Vec<u8>) だと他の箇所ももっと楽にできたりしないかな
    // 復元する時にbit化前行列が単位行列の一行だとそのままコピーで済んだりとか
    // 今はcopyが1/8 * 8になってて、直接1やった方が速い気もする
    // schedulerのoptimizedでなんとかならないもんかねえ
    // schedulerにはmatrixそのまま渡してるわけだから……
    for i in 0..parity.len() {
        // mslice[width * i..width * (i + 1)].copy_from_slice(&parity[i]);
        fast_array_copy(&mut mslice[width * i..width * (i + 1)], &parity[i]);
    }
    for i in 0..original_height {
        if i != last_idx_of_erased_data {
            let slice = &output[original_width * i..original_width * (i + 1)];
            fast_array_xor(mslice, slice);
        }
    }

    output
}

pub fn to_bitmatrix4(
    commands: Vec<Command>,
    nr_data: usize,
    nr_parity: usize,
    data: &[&[u8]],
) -> Vec<Vec<u8>> {
    let mut output: Vec<Vec<u8>> = Vec::new();
    let width = data[0].len();

    for _ in 0..(nr_data + nr_parity) * 8 {
        let mut v = Vec::with_capacity(width);
        unsafe {
            v.set_len(width);
        }
        output.push(v);
    }

    // コマンドを解釈していく
    for c in commands {
        if let Command::Copy(from, to) = c {
            output[to].copy_from_slice(data[from]);
        } else if let Command::Xor(from, to) = c {
            xor_vec(&mut output[to], data[from]);
        }
    }

    output
}

pub fn to_bitmatrix4_optim(
    commands: Vec<Command>,
    nr_data: usize,
    nr_parity: usize,
    data: &[&[u8]],
) -> Vec<Vec<u8>> {
    let mut output: Vec<Vec<u8>> = Vec::new();
    let width = data[0].len();

    for _ in 0..(nr_data + nr_parity) * 8 {
        let mut v = Vec::with_capacity(width);
        unsafe {
            v.set_len(width);
        }
        output.push(v);
    }

    // コマンドを解釈していく
    for c in commands {
        if let Command::Copy(from, to) = c {
            fast_array_copy(&mut output[to], data[from]);
        } else if let Command::Xor(from, to) = c {
            fast_array_xor(&mut output[to], data[from]);
        }
    }

    output
}

pub fn to_bitmatrix3(m: &Matrix<GF_2_8>, data: &[&[u8]]) -> Vec<Vec<u8>> {
    // まず (m.height() * 8, data.width()) 行列を作る
    let mut output: Vec<Vec<u8>> = Vec::new();
    let width = data[0].len();
    for _ in 0..m.height() * 8 {
        let mut v = Vec::with_capacity(width);
        unsafe {
            v.set_len(width);
        }
        output.push(v);
    }

    let commands = matrix_to_commands(m);

    // コマンドを解釈していく
    for c in commands {
        if let Command::Copy(from, to) = c {
            output[to].copy_from_slice(data[from]);
        } else if let Command::Xor(from, to) = c {
            xor_vec(&mut output[to], data[from]);
        }
    }

    output
}

pub fn to_bitmatrix2(m: &Matrix<GF_2_8>, data: &[&[u8]]) -> Vec<Vec<u8>> {
    let width = data[0].len();

    let mut output: Vec<Vec<u8>> = Vec::new();

    for i in 0..m.height() {
        let row: Vec<u8> = m[i].as_vec().iter().map(|e| u8::from(*e)).collect();
        let v: Vec<[u8; 8]> = row.iter().map(|e| BIT_GF_2_8_IMPL.fast_conv(*e)).collect();
        let v: Vec<Vec<u8>> = expand(&v);

        for (j, e) in v.iter().enumerate() {
            let mut rowv: Vec<u8> = Vec::with_capacity(width);
            unsafe {
                rowv.set_len(width);
            }
            let commands = gen_command(e, i * 8 + j);
            for c in commands {
                if let Command::Copy(from, _) = c {
                    rowv.copy_from_slice(data[from]);
                } else if let Command::Xor(from, _) = c {
                    xor_vec(&mut rowv, data[from]);
                }
            }
            output.push(rowv);
        }
    }

    output
}

pub fn to_bitmatrix(m: &Matrix<GF_2_8>, data: &[&[u8]]) -> Vec<Vec<u8>> {
    let mut output = Vec::new();

    for i in 0..m.height() {
        // 行ベクトルを取り出す
        let row: Vec<u8> = m[i].as_vec().iter().map(|e| u8::from(*e)).collect();
        // 縦にふくらませる
        let v: Vec<[u8; 8]> = row.iter().map(|e| BIT_GF_2_8_IMPL.fast_conv(*e)).collect();
        let v: Vec<Vec<u8>> = expand(&v);

        for e in v {
            output.push(bit_matrix_row_prod(&e, &data));
        }
    }

    output
}

pub fn encode_by_bit(m: &Matrix<GF_2_8>, data: &[u8]) -> Vec<Vec<u8>> {
    let data = vec_to_quasi_matrix(&data, m.width() * 8);

    to_bitmatrix(&m, &data)
}

#[allow(non_camel_case_types)]
pub struct Bit_GF_2_8_impl {
    ppoly: Poly<GF_2>,
    table: Vec<[u8; 8]>,
}

fn xors(x: u8) -> u8 {
    (x.count_ones() % 2) as u8
}

impl Bit_GF_2_8_impl {
    pub fn ppoly(&self) -> Poly<GF_2> {
        self.ppoly.clone()
    }

    pub fn new(ppoly: Poly<GF_2>) -> Self {
        Self {
            ppoly,
            table: Vec::new(),
        }
    }

    pub fn fast_conv(&self, u: u8) -> [u8; 8] {
        self.table[u as usize]
    }

    pub fn setup(&mut self) {
        for i in 0u8..=255u8 {
            let t = self.conv(i.into());
            self.table.push(t);
        }
    }

    pub fn build(ppoly: Poly<GF_2>) -> Self {
        let mut imp = Self::new(ppoly);
        imp.setup();
        imp
    }

    pub fn mul(&self, p: GF_2_8, q: GF_2_8) -> GF_2_8 {
        /*
         * p[0] = c_of_deg7
         * p[1] = c_of_deg6
         * ..
         * p[7] = c_of_deg0
         * の並びに注意
         *
         * (p[0]) q[deg7]
         * (p[1]) q[deg6]
         * (....) ...
         * (p[7]) q[deg0]
         * で行列積する感じ
         */

        let p: [u8; 8] = self.table[u8::from(p) as usize];

        // (q >> deg) = 係数
        let q: u8 = q.into();

        let mut r: u8 = 0;

        for (i, pe) in p.iter().enumerate() {
            let v: u8 = pe & q;
            let v: u8 = xors(v);

            let deg = 7 - i;
            r |= v << deg;
        }
        r.into()
    }

    // matrix積の布石
    pub fn mul2(&self, p: GF_2_8, q: GF_2_8) -> GF_2_8 {
        let p: [u8; 8] = self.table[u8::from(p) as usize];

        let mut r: u8 = 0;

        for (i, row) in p.iter().enumerate() {
            let deg1 = 7 - i;

            for j in 0..8 {
                let deg2 = 7 - j;
                if (row >> deg2) & 1 == 1 {
                    r ^= q.coef(deg2) << deg1;
                } else {
                    // 0かけてxorなので何もしなくて良い
                }
            }
        }

        r.into()
    }

    /*
     * p = c_7 x^7 + c_6 x^6 + ... + c_1 x^1 + c_0
     * について
     * (p q) % ppoly (= x^8 + x^4 + x^3 + x^2 + 1)
     * を行列演算で計算するための部分最適化をする
     *
     * (p q) % ppoly の x^7 の係数は
     *
     * (p * d_7 x^7) % ppoly@x^7 + (p * d_6 x^6) % ppoly@x^7 + ... + (p * d_0) * ppoly@x^7
     *
     * (p * d_7 x^7) % ppol = d_7 (p * x^7) % ppoly
     * なので (p * x^7) % ppoly を先に計算しておけば良い
     */
    pub fn conv(&self, p: GF_2_8) -> [u8; 8] {
        // r[degree]
        // この並びにしておくと conv(1) が単位行列になって嬉しい
        let mut r: [u8; 8] = [7, 6, 5, 4, 3, 2, 1, 0];

        // 次数iの計算
        // assume: rv = &r[i]
        for i in 0..8 {
            let i: u32 = i as u32;
            let mut v: u8 = 0;
            // 次数deg毎に積 p をとり
            // 次数iの係数を取り出しておく
            for deg in 0..8 {
                let p_ = (p.to_poly() * Poly::from_mono(deg, GF_2::ONE)) % self.ppoly.clone();

                /*
                println!("{} * {} % {} = {}",
                         p.to_poly().to_string_as_poly(),
                         Poly::from_mono(deg, GF_2::ONE).to_string_as_poly(),
                         self.ppoly.clone().to_string_as_poly(),
                         p_.to_string_as_poly());
                 */

                v |= p_.at(&i).to_u8() << deg;
            }
            r[(7 - i) as usize] = v;
        }

        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vandermonde::*;
    use std::time::Instant;

    fn ppoly() -> Poly<GF_2> {
        Poly::from_vec(vec![
            (8, GF_2::ONE),
            (4, GF_2::ONE),
            (3, GF_2::ONE),
            (2, GF_2::ONE),
            (0, GF_2::ONE),
        ])
    }

    #[test]
    fn test_gen_command() {
        use super::Command::*;

        let to = 0;
        let v = vec![0b10000000];
        let c = gen_command(&v, to);

        assert!(c == vec![Copy(0, to)]);

        let v = vec![0b00100000];
        let c = gen_command(&v, to);

        assert!(c == vec![Copy(2, to)]);

        let to = 5;
        let v = vec![0b10100000];
        let c = gen_command(&v, to);

        assert!(c == vec![Copy(0, to), Xor(2, to)]);

        let to = 5;
        let v = vec![0b10100000, 0b01000000];
        let c = gen_command(&v, to);

        assert!(c == vec![Copy(0, to), Xor(2, to), Xor(9, to)]);

        let to = 5;
        let v = vec![0b10100000, 0b01000001];
        let c = gen_command(&v, to);

        assert!(c == vec![Copy(0, to), Xor(2, to), Xor(9, to), Xor(15, to)]);

        let to = 42;
        let v = vec![0b00000000, 0b01000000];
        let c = gen_command(&v, to);

        assert!(c == vec![Copy(9, to)]);
    }

    #[test]
    fn test_expand() {
        let v = vec![[0, 1, 2, 3, 4, 5, 6, 7], [10, 11, 12, 13, 14, 15, 16, 17]];

        let w = vec![
            vec![0, 10],
            vec![1, 11],
            vec![2, 12],
            vec![3, 13],
            vec![4, 14],
            vec![5, 15],
            vec![6, 16],
            vec![7, 17],
        ];

        assert!(expand(&v) == w);
    }

    #[test]
    fn test_make_sched() {
        let v = vec![0b00000101, 0b00000010];

        assert!(make_sched(&v) == vec![5, 7, 14]);
    }

    #[test]
    fn bit_mul_test1() {
        for i in 1..=100 {
            let idm: Matrix<GF_2_8> = Matrix::identity(i);
            let data: Vec<u8> = (0..8 * i).map(|_| rand::random::<u8>()).collect();

            let enc_data: Vec<Vec<u8>> = encode_by_bit(&idm, &data);
            let dec_data = enc_data.concat();
            assert!(data == dec_data);
        }
    }

    #[test]
    fn bit_mul_test2() {
        for i in 1..=100 {
            let idm: Matrix<GF_2_8> = Matrix::identity(i);
            let data: Vec<u8> = (0..100 * (8 * i)).map(|_| rand::random::<u8>()).collect();

            let enc_data: Vec<Vec<u8>> = encode_by_bit(&idm, &data);
            let dec_data = enc_data.concat();
            assert!(data == dec_data);
        }
    }

    #[test]
    fn bit_mul_test3() {
        for i in 1..=100 {
            let velems: Vec<GF_2_8> = (1..=i)
                .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
                .collect();

            let mut m: Matrix<GF_2_8> = vandermonde(
                MatrixSize {
                    height: i,
                    width: i,
                },
                &velems,
            )
            .unwrap();

            let data: Vec<u8> = (0..100 * (8 * i)).map(|_| rand::random::<u8>()).collect();

            let enc_data: Vec<Vec<u8>> = encode_by_bit(&m, &data);

            let inv = m.inverse().unwrap();

            let enc_data: Vec<&[u8]> = enc_data.iter().map(|v| &v[..]).collect();
            let dec_data: Vec<Vec<u8>> = to_bitmatrix(&inv, &enc_data);
            let dec_data = dec_data.concat();

            assert!(data == dec_data);
        }
    }

    #[test]
    fn enc_dec_test1() {
        let nr_data = 2;
        let nr_parity = 1;

        let velems: Vec<GF_2_8> = (1..=nr_data + nr_parity)
            .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
            .collect();

        let mut m: Matrix<GF_2_8> = vandermonde(
            MatrixSize {
                height: nr_data + nr_parity,
                width: nr_data,
            },
            &velems,
        )
        .unwrap();

        let data: Vec<u8> = (0..8 * nr_data).map(|_| rand::random::<u8>()).collect();

        let mut enc_data: Vec<Vec<u8>> = {
            encode_by_bit(&m, &data)

            // memory_optimized_mul(&m, &vec_to_quasi_matrix(&data, nr_data))
        };

        println!(
            "enc height = {}, width = {}\n",
            enc_data.len(),
            enc_data[0].len()
        );

        // データブロックを復旧するパターン
        let remove_col = 2;
        for i in (0..8).rev() {
            enc_data.remove(remove_col * 8 + i);
        }

        m.drop_columns(vec![remove_col]);
        let inv = m.inverse().unwrap();

        let enc_data: Vec<&[u8]> = enc_data.iter().map(|v| &v[..]).collect();

        let dec_data: Vec<Vec<u8>> = {
            to_bitmatrix(&inv, &enc_data)

            // memory_optimized_mul(&inv, &enc_data)
        };

        let dec_data = dec_data.concat();

        assert!(data == dec_data);
    }

    #[test]
    fn test_to_bitmatrix2() {
        let nr_data = 20;
        let nr_parity = 10;

        let velems: Vec<GF_2_8> = (1..=nr_data + nr_parity)
            .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
            .collect();

        let m: Matrix<GF_2_8> = vandermonde(
            MatrixSize {
                height: nr_data + nr_parity,
                width: nr_data,
            },
            &velems,
        )
        .unwrap();

        let data: Vec<u8> = (0..8 * nr_data).map(|_| rand::random::<u8>()).collect();

        let enc_data1: Vec<Vec<u8>> = to_bitmatrix(&m, &vec_to_quasi_matrix(&data, 8 * nr_data));
        let enc_data2: Vec<Vec<u8>> = to_bitmatrix2(&m, &vec_to_quasi_matrix(&data, 8 * nr_data));

        assert!(enc_data1 == enc_data2);
    }

    #[test]
    fn test_to_bitmatrix3() {
        let nr_data = 20;
        let nr_parity = 10;

        let velems: Vec<GF_2_8> = (1..=nr_data + nr_parity)
            .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
            .collect();

        let m: Matrix<GF_2_8> = vandermonde(
            MatrixSize {
                height: nr_data + nr_parity,
                width: nr_data,
            },
            &velems,
        )
        .unwrap();

        let data: Vec<u8> = (0..8 * nr_data).map(|_| rand::random::<u8>()).collect();

        let enc_data1: Vec<Vec<u8>> = to_bitmatrix(&m, &vec_to_quasi_matrix(&data, 8 * nr_data));
        let enc_data3: Vec<Vec<u8>> = to_bitmatrix3(&m, &vec_to_quasi_matrix(&data, 8 * nr_data));

        assert!(enc_data1 == enc_data3);
    }

    #[test]
    fn test_bitenc_dec1() {
        let nr_data = 10;
        let nr_parity = 4;

        let velems: Vec<GF_2_8> = (1..=nr_data + nr_parity)
            .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
            .collect();

        let m: Matrix<GF_2_8> = modified_systematic_vandermonde(
            MatrixSize {
                height: nr_data + nr_parity,
                width: nr_data,
            },
            &velems,
        )
        .unwrap();

        let data: Vec<u8> = (0..8 * nr_data).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc1(&m, &data);

        encoded.remove(3);
        encoded.remove(2);
        encoded.remove(1);
        encoded.remove(0);

        let decoded = bitmatrix_dec1(&m, &encoded);

        assert!(data == decoded);
    }

    #[test]
    fn mul_test() {
        let i1 = GF_2_8_impl::new(ppoly());
        let i2 = Bit_GF_2_8_impl::build(ppoly());

        for i in 0u8..=255u8 {
            for j in 0u8..=255u8 {
                let r1 = i1.mul(i.into(), j.into());
                let r2 = i2.mul(i.into(), j.into());
                let r3 = i2.mul2(i.into(), j.into());

                assert!(r1 == r2, "[1] {} * {}", i, j);
                assert!(r2 == r3, "[2] {} * {}", i, j);
            }
        }
    }

    #[test]
    fn speed_test() {
        let i1 = GF_2_8_impl::new(ppoly());
        let i2 = Bit_GF_2_8_impl::build(ppoly());

        let mut r1: u8 = 0;
        let mut r2: u8 = 0;
        let mut r3: u8 = 0;

        let t1 = Instant::now();

        for _ in 0..100 {
            for i in 0u8..=255u8 {
                for j in 0u8..=255u8 {
                    r1 ^= u8::from(i1.mul(i.into(), j.into()));
                }
            }
        }
        println!("naiive mul = {:?}", t1.elapsed());

        let t2 = Instant::now();
        for _ in 0..100 {
            for i in 0u8..=255u8 {
                for j in 0u8..=255u8 {
                    r2 ^= u8::from(i2.mul(i.into(), j.into()));
                }
            }
        }
        println!("bitmatrix mul = {:?}", t2.elapsed());

        let t3 = Instant::now();
        for _ in 0..100 {
            for i in 0u8..=255u8 {
                for j in 0u8..=255u8 {
                    r3 ^= u8::from(i2.mul2(i.into(), j.into()));
                }
            }
        }
        println!("bitmatrix mul2 = {:?}", t3.elapsed());

        assert!(r1 == r2);
        assert!(r2 == r3);
    }
}
