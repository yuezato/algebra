use crate::bit_based_gf8::*;
use crate::erasure_code::*;
use crate::fin_field::*;
use crate::matrix::*;
use fast_array_ops::*;
use std::convert::TryInto;

/*
 * Paddingは工夫しないとメモリコピーを発生させてしまう。
 *
 * 仮想的に 0 埋めした行ベクトルを表現するための構造体
 * Rowを利用する。
 *
 * x x x x x x x
 * x x x x x x x
 * x x x . . . .
 * . . . . . . .
 * x x 0 0 0 0 0
 * 0 0 0 0 0 0 0
 * 0 0 0 0 0 0 0
 */

#[derive(Debug)]
pub struct Row<'a> {
    len: usize, // 有効なエントリ数
    inner: &'a [u8],
}

pub fn split_and_padding<'a>(v: &'a [u8], height: usize) -> Vec<Row<'a>> {
    let mut r = Vec::new();
    let mut row_elems = v.len() / height;

    if v.len() % height != 0 {
        row_elems += 1; // 割り切れないのでpaddingする
    }

    let mut residual = v.len();
    let mut cur = 0;

    for _ in 0..height {
        if residual >= row_elems {
            r.push(Row {
                len: row_elems,
                inner: &v[cur..cur + row_elems],
            });
            residual -= row_elems;
            cur += row_elems;
        } else {
            r.push(Row {
                len: residual,
                inner: &v[cur..cur + residual],
            });
            residual = 0;
        }
    }

    r
}

pub fn normal_array_copy(dst: &mut [u8], src: &[u8]) {
    dst.copy_from_slice(src);
}

pub fn is_aligned<T>(p: *const T, a: usize) -> bool {
    (p as usize) % a == 0
}

// alignmentも長さも何も期待できない場合の実装
pub fn normal_array_xor(dst: &mut [u8], src: &[u8]) {
    if is_aligned(dst.as_ptr(), 32) && is_aligned(src.as_ptr(), 32) && dst.len() % 64 == 0 {
        fast_array_xor(dst, src);
    } else {
        for i in 0..src.len() {
            dst[i] ^= src[i];
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct EncInfo {
    block_id: u8,
    original_len: u32,
    nr_data: u8,
    nr_parity: u8,
}

fn cast(slice: &[u8]) -> &[u8; 7] {
    let ptr = slice.as_ptr() as *const [u8; 7];
    unsafe { &*ptr }
}

impl EncInfo {
    pub fn to_bytes(&self) -> [u8; 7] {
        let mut v = [0, 0, 0, 0, 0, 0, 0];
        v[0] = self.block_id;
        v[1] = (self.original_len & 0xff) as u8;
        v[2] = ((self.original_len >> 8) & 0xff) as u8;
        v[3] = ((self.original_len >> 16) & 0xff) as u8;
        v[4] = ((self.original_len >> 24) & 0xff) as u8;
        v[5] = self.nr_data;
        v[6] = self.nr_parity;
        v
    }

    pub fn from_bytes(bytes: &[u8; 7]) -> EncInfo {
        EncInfo {
            block_id: bytes[0],
            original_len: u32::from_be_bytes([bytes[4], bytes[3], bytes[2], bytes[1]]),
            nr_data: bytes[5],
            nr_parity: bytes[6],
        }
    }
}

pub fn to_bitmatrix_pad<'a>(
    commands: Vec<Command>,
    nr_data: usize,
    nr_parity: usize,
    data: &[Row<'a>],
    original_len: u32,
) -> Vec<Vec<u8>> {
    let mut output: Vec<Vec<u8>> = Vec::new();
    let width = data[0].len;

    // 符号化済みデータ行列の初期化
    for block_id in 0..(nr_data + nr_parity) {
        // zeroing is important wrt padding
        // 零化しないと失敗する理屈も書いておきたい
        // でも本来零であるところに謎の値があって
        // それが他の値を計算するところにも使われるんだから
        // 当たり前といえば当たり前か
        let mut v: Vec<u8> = vec![0; width * 8 + 7];
        let bytes = (EncInfo {
            nr_data: nr_data.try_into().unwrap(),
            nr_parity: nr_parity.try_into().unwrap(),
            original_len,
            block_id: block_id.try_into().unwrap(),
        })
        .to_bytes();
        v[width * 8 + 0] = bytes[0];
        v[width * 8 + 1] = bytes[1];
        v[width * 8 + 2] = bytes[2];
        v[width * 8 + 3] = bytes[3];
        v[width * 8 + 4] = bytes[4];
        v[width * 8 + 5] = bytes[5];
        v[width * 8 + 6] = bytes[6];

        output.push(v);
    }

    // コマンドを解釈していく
    for c in commands {
        if let Command::Copy(from, to) = c {
            let len = data[from].len;

            let src = &data[from].inner[0..len];

            let to1 = to / 8;
            let to2 = to % 8;
            let dst = &mut output[to1][to2 * width..to2 * width + len];

            normal_array_copy(dst, src);
        } else if let Command::Xor(from, to) = c {
            let len = data[from].len;

            let src = &data[from].inner[0..len];

            let to1 = to / 8;
            let to2 = to % 8;
            let dst = &mut output[to1][to2 * width..to2 * width + len];

            normal_array_xor(dst, src);
        }
    }

    output
}

pub fn to_bitmatrix_pad_optim<'a>(
    commands: Vec<Command>,
    nr_data: usize,
    nr_parity: usize,
    data: &[Row<'a>],
    orig_data: &[u8],
    original_len: u32,
) -> Vec<Vec<u8>> {
    let mut output: Vec<Vec<u8>> = Vec::new();
    let width = data[0].len;

    // 符号化済みデータ行列の初期化
    let mut cur = 0;
    let mut residual = original_len as usize;
    for block_id in 0..nr_data {
        let mut v: Vec<u8> = vec![0; width * 8 + 7];       

        if residual >= width * 8 {
            normal_array_copy(&mut v[0..width*8], &orig_data[cur..cur+width*8]);
            residual -= width * 8;
            cur += width * 8;
        } else {
            normal_array_copy(&mut v[0..residual], &orig_data[cur..cur+residual]);
            residual = 0;
        }
        
        let bytes = (EncInfo {
            nr_data: nr_data.try_into().unwrap(),
            nr_parity: nr_parity.try_into().unwrap(),
            original_len,
            block_id: block_id.try_into().unwrap(),
        })
        .to_bytes();
        v[width * 8 + 0] = bytes[0];
        v[width * 8 + 1] = bytes[1];
        v[width * 8 + 2] = bytes[2];
        v[width * 8 + 3] = bytes[3];
        v[width * 8 + 4] = bytes[4];
        v[width * 8 + 5] = bytes[5];
        v[width * 8 + 6] = bytes[6];

        output.push(v);
    }
    
    for block_id in nr_data..nr_data+nr_parity {
        let mut v: Vec<u8> = vec![0; width * 8 + 7];
        let bytes = (EncInfo {
            nr_data: nr_data.try_into().unwrap(),
            nr_parity: nr_parity.try_into().unwrap(),
            original_len,
            block_id: block_id.try_into().unwrap(),
        })
        .to_bytes();
        v[width * 8 + 0] = bytes[0];
        v[width * 8 + 1] = bytes[1];
        v[width * 8 + 2] = bytes[2];
        v[width * 8 + 3] = bytes[3];
        v[width * 8 + 4] = bytes[4];
        v[width * 8 + 5] = bytes[5];
        v[width * 8 + 6] = bytes[6];

        output.push(v);
    }

    // コマンドを解釈していく
    // ただし to が data の場合には既に完了しているので
    // 何もしない
    for c in commands {
        if c.to() / 8 < nr_data {
            continue;
        }
        if let Command::Copy(from, to) = c {
            let len = data[from].len;

            let src = &data[from].inner[0..len];

            let to1 = to / 8;
            let to2 = to % 8;
            let dst = &mut output[to1][to2 * width..to2 * width + len];

            normal_array_copy(dst, src);
        } else if let Command::Xor(from, to) = c {
            let len = data[from].len;

            let src = &data[from].inner[0..len];

            let to1 = to / 8;
            let to2 = to % 8;
            let dst = &mut output[to1][to2 * width..to2 * width + len];

            normal_array_xor(dst, src);
        }
    }

    output
}

pub fn bitmatrix_enc_with_info(
    m: &Matrix<GF_2_8>,
    nr_data: usize,
    nr_parity: usize,
    data: &[u8],
) -> Vec<Vec<u8>> {
    let commands = matrix_to_commands(&m);
    let splitted = split_and_padding(&data, 8 * nr_data);

    to_bitmatrix_pad(
        commands,
        nr_data,
        nr_parity,
        &splitted,
        data.len().try_into().unwrap(),
    )
}

pub fn bitmatrix_enc_with_info_optim(
    m: &Matrix<GF_2_8>,
    nr_data: usize,
    nr_parity: usize,
    data: &[u8],
) -> Vec<Vec<u8>> {
    let commands = matrix_to_commands(&m);
    let splitted = split_and_padding(&data, 8 * nr_data);

    to_bitmatrix_pad_optim(
        commands,
        nr_data,
        nr_parity,
        &splitted,
        data,
        data.len().try_into().unwrap(),
    )
}

pub fn retrieve_info(data: &Vec<Vec<u8>>) -> (EncInfo, Vec<usize>) {
    let l = data[0].len();
    let e = EncInfo::from_bytes(cast(&data[0][l - 7..l]));

    let mut alive = Vec::new();

    for d in data {
        let e_ = EncInfo::from_bytes(cast(&d[l - 7..l]));
        let block = e_.block_id;
        alive.push(block as usize);
    }

    (e, alive)
}

pub fn to_bitmatrix_pad2(commands: Vec<Command>, data: &[&[u8]]) -> Vec<u8> {
    // for i in 0..data.len() { dbg!((i, data[i].len())); }

    let width = data[0].len();
    let orig_width = width / 8;
    let padded_len = width * data.len();

    let mut output: Vec<u8> = Vec::with_capacity(padded_len);
    unsafe {
        output.set_len(padded_len);
    }

    // コマンドを解釈していく
    for c in commands {
        if let Command::Copy(from, to) = c {
            let from1 = from / 8;
            let from2 = from % 8;
            let src = &data[from1][from2 * orig_width..(from2 + 1) * orig_width];

            let dst = &mut output[to * orig_width..(to + 1) * orig_width];
            normal_array_copy(dst, src);
        } else if let Command::Xor(from, to) = c {
            let from1 = from / 8;
            let from2 = from % 8;
            let src = &data[from1][from2 * orig_width..(from2 + 1) * orig_width];

            let dst = &mut output[to * orig_width..(to + 1) * orig_width];
            normal_array_xor(dst, src);
        }
    }
    output
}

pub fn bitmatrix_dec_with_info(m: &Matrix<GF_2_8>, data: &Vec<Vec<u8>>) -> Vec<u8> {
    let (enc_info, alive_ids) = retrieve_info(&data);

    let alive =
        AliveBlocks::from_alive_vec((enc_info.nr_data + enc_info.nr_parity).into(), &alive_ids);

    let invm = decode_matrix(m.clone(), &alive).unwrap();

    let commands = matrix_to_commands(&invm);

    let mut data2: Vec<&[u8]> = Vec::new();
    for v in data {
        data2.push(&v[0..v.len() - 7]);
    }

    let mut v = to_bitmatrix_pad2(commands, &data2);

    unsafe {
        v.set_len(enc_info.original_len.try_into().unwrap());
    }

    v
}

pub fn bitmatrix_dec_with_info_optim(m: &Matrix<GF_2_8>, data: &Vec<Vec<u8>>) -> Vec<u8> {
    let (enc_info, alive_ids) = retrieve_info(&data);

    let alive =
        AliveBlocks::from_alive_vec((enc_info.nr_data + enc_info.nr_parity).into(), &alive_ids);

    let invm = decode_matrix(m.clone(), &alive).unwrap();

    let commands = matrix_to_commands(&invm);

    let mut data2: Vec<&[u8]> = Vec::new();
    let mut parity = None;
    for v in data {
        data2.push(&v[0..v.len() - 7]);
        let e = EncInfo::from_bytes(cast(&v[v.len() - 7..v.len()]));
        if e.block_id == e.nr_data {
            parity = Some(&v[0..v.len() - 7]);
        }
    }

    let erased_data_ids: Vec<u8> = (0..enc_info.nr_data)
        .filter(|i| {
            let i = *i as usize;
            !alive_ids.contains(&i)
        })
        .collect();

    let mut v = if erased_data_ids.is_empty() {
        // データブロックが無傷の場合
        let width = data2[0].len();
        let size = data2.len() * data2[0].len();
        let mut v = Vec::with_capacity(size);
        unsafe {
            v.set_len(size);
        }

        for i in 0..enc_info.nr_data {
            let i = i as usize;
            v[i * width..(i + 1) * width].copy_from_slice(data2[i]);
        }

        v
    } else if let Some(parity) = parity {
        // topmost parityが生き残っている場合
        let an_erased_data_id: usize = *erased_data_ids.first().unwrap() as usize;

        let mut filtered_commands = Vec::new();
        for c in commands {
            if !(an_erased_data_id * 8..(an_erased_data_id + 1) * 8).contains(&c.to()) {
                filtered_commands.push(c);
            } else {
                println!("{:?}", c);
            }
        }
        with_topmost_parity(filtered_commands, &data2, &parity, an_erased_data_id)
    } else {
        // topmost parityがない通常の場合
        to_bitmatrix_pad2(commands, &data2)
    };

    unsafe {
        v.set_len(enc_info.original_len.try_into().unwrap());
    }

    v
}

pub fn with_topmost_parity(
    filtered_commands: Vec<Command>,
    data: &[&[u8]], // embedded infoは除去済み
    parity: &[u8],  // embedded infoは除去済み
    last_idx_of_erased_data: usize,
) -> Vec<u8> {
    let width = data[0].len();
    let orig_width = width / 8;
    let height = data.len();

    let mut output: Vec<u8> = Vec::with_capacity(width * height);
    unsafe { output.set_len(width * height) }

    for c in filtered_commands {
        if let Command::Copy(from, to) = c {
            let from1 = from / 8;
            let from2 = from % 8;
            let src = &data[from1][from2 * orig_width..(from2 + 1) * orig_width];

            let dst = &mut output[to * orig_width..(to + 1) * orig_width];
            normal_array_copy(dst, src);
        } else if let Command::Xor(from, to) = c {
            let from1 = from / 8;
            let from2 = from % 8;
            let src = &data[from1][from2 * orig_width..(from2 + 1) * orig_width];

            let dst = &mut output[to * orig_width..(to + 1) * orig_width];
            normal_array_xor(dst, src);
        }
    }

    let mslice: &mut [u8] = unsafe {
        let ptr: *mut u8 = output.as_ptr() as *mut u8;
        let ptr: *mut u8 = ptr.add(width * last_idx_of_erased_data);
        std::slice::from_raw_parts_mut(ptr, width)
    };

    normal_array_copy(mslice, parity);

    for i in 0..height {
        if i != last_idx_of_erased_data {
            let slice = &output[width * i..width * (i + 1)];
            normal_array_xor(mslice, slice);
        }
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::*;
    use crate::vandermonde::*;

    #[test]
    fn enc_dec_test1() {
        let idm: Matrix<GF_2_8> = Matrix::identity(3);
        let data: Vec<u8> = (0..100).map(|_| rand::random::<u8>()).collect();

        let encoded = bitmatrix_enc_with_info(&idm, 3, 0, &data);

        dbg!(encoded[0].len());
        dbg!(encoded[1].len());
        dbg!(encoded[2].len());

        let decoded = bitmatrix_dec_with_info(&idm, &encoded);

        dbg!(&data);
        dbg!(&decoded);

        assert!(data == decoded);
    }

    #[test]
    fn enc_dec_test2() {
        let nr_data = 3;
        let nr_parity = 1;

        let velems: Vec<GF_2_8> = (1..=nr_data + nr_parity)
            .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
            .collect();

        let mut vanderm: Matrix<GF_2_8> = modified_systematic_vandermonde(
            MatrixSize {
                height: nr_data + nr_parity,
                width: nr_data,
            },
            &velems,
        )
        .unwrap();

        let data: Vec<u8> = (0..48).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

        encoded.remove(0);

        let decoded = bitmatrix_dec_with_info(&vanderm, &encoded);

        dbg!(&data);
        dbg!(&decoded);

        assert!(data == decoded);
    }

    fn make_vandermonde_matrix(nr_data: usize, nr_parity: usize) -> Matrix<GF_2_8> {
        let velems: Vec<GF_2_8> = (1..=nr_data + nr_parity)
            .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
            .collect();

        modified_systematic_vandermonde(
            MatrixSize {
                height: nr_data + nr_parity,
                width: nr_data,
            },
            &velems,
        )
        .unwrap()
    }

    #[test]
    fn enc_dec_test3() {
        let nr_data = 10;
        let nr_parity = 4;

        let velems: Vec<GF_2_8> = (1..=nr_data + nr_parity)
            .map(|i| GF_2_8::PRIMITIVE_ELEMENT.exp(i as u32))
            .collect();

        let mut vanderm: Matrix<GF_2_8> = modified_systematic_vandermonde(
            MatrixSize {
                height: nr_data + nr_parity,
                width: nr_data,
            },
            &velems,
        )
        .unwrap();

        let data: Vec<u8> = (0..93 * 1024).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

        encoded.remove(12);
        encoded.remove(10);
        encoded.remove(3);
        encoded.remove(0);

        let decoded = bitmatrix_dec_with_info(&vanderm, &encoded);

        assert!(data == decoded);
    }

    #[test]
    fn enc_optimdec_test1() {
        let nr_data = 10;
        let nr_parity = 4;

        let vanderm = make_vandermonde_matrix(nr_data, nr_parity);

        let data: Vec<u8> = (0..93 * 1024).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

        // データが無傷のパターン
        encoded.remove(13);
        encoded.remove(12);
        encoded.remove(11);
        encoded.remove(10);

        let decoded = bitmatrix_dec_with_info_optim(&vanderm, &encoded);

        assert!(data == decoded);
    }

    #[test]
    fn enc_optimdec_test2() {
        let nr_data = 10;
        let nr_parity = 4;

        let vanderm = make_vandermonde_matrix(nr_data, nr_parity);

        let data: Vec<u8> = (0..93 * 1024).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

        // データもtopmost parityも消失している場合
        encoded.remove(13);
        encoded.remove(12);
        encoded.remove(10);
        encoded.remove(9);

        let decoded = bitmatrix_dec_with_info_optim(&vanderm, &encoded);

        assert!(data == decoded);
    }

    #[test]
    fn enc_optimdec_test3() {
        let nr_data = 10;
        let nr_parity = 4;

        let vanderm = make_vandermonde_matrix(nr_data, nr_parity);

        let data: Vec<u8> = (0..93 * 1024).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

        // データが1つ消失
        // topmost parityは生き残っている
        encoded.remove(13);
        encoded.remove(12);
        encoded.remove(11);
        encoded.remove(9);

        let decoded = bitmatrix_dec_with_info_optim(&vanderm, &encoded);

        assert!(data == decoded);
    }

    #[test]
    fn enc_optimdec_test4() {
        let nr_data = 10;
        let nr_parity = 4;

        let vanderm = make_vandermonde_matrix(nr_data, nr_parity);

        let data: Vec<u8> = (0..93 * 1024).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

        // データが2つ消失
        // topmost parityは生き残っている
        encoded.remove(13);
        encoded.remove(12);
        encoded.remove(9);
        encoded.remove(8);

        let decoded = bitmatrix_dec_with_info_optim(&vanderm, &encoded);

        assert!(data == decoded);
    }

    #[test]
    fn enc_optimdec_test5() {
        let nr_data = 10;
        let nr_parity = 4;

        let vanderm = make_vandermonde_matrix(nr_data, nr_parity);

        let data: Vec<u8> = (0..93 * 1024).map(|_| rand::random::<u8>()).collect();

        let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

        // データが3つ消失
        // topmost parityは生き残っている
        encoded.remove(13);
        encoded.remove(9);
        encoded.remove(8);
        encoded.remove(7);

        let decoded = bitmatrix_dec_with_info_optim(&vanderm, &encoded);

        assert!(data == decoded);
    }

    #[test]
    fn enc_optimdec_test6() {
        let nr_data = 10;
        let nr_parity = 4;

        let vanderm = make_vandermonde_matrix(nr_data, nr_parity);

        for i in 1..100 {
            let data: Vec<u8> = (0..i * 1024).map(|_| rand::random::<u8>()).collect();

            let mut encoded = bitmatrix_enc_with_info(&vanderm, nr_data, nr_parity, &data);

            // データが4つ消失
            encoded.remove(3);
            encoded.remove(2);
            encoded.remove(1);
            encoded.remove(0);

            let decoded = bitmatrix_dec_with_info_optim(&vanderm, &encoded);

            assert_eq!(data, decoded);
        }
    }

    #[test]
    fn optimenc_optimdec_test1() {
        let nr_data = 10;
        let nr_parity = 4;

        let vanderm = make_vandermonde_matrix(nr_data, nr_parity);

        for i in 1..100 {
            let data: Vec<u8> = (0..i * 1024).map(|_| rand::random::<u8>()).collect();

            let mut encoded = bitmatrix_enc_with_info_optim(&vanderm, nr_data, nr_parity, &data);

            // データが4つ消失
            encoded.remove(3);
            encoded.remove(2);
            encoded.remove(1);
            encoded.remove(0);

            let decoded = bitmatrix_dec_with_info_optim(&vanderm, &encoded);

            assert_eq!(data, decoded);
        }
    }

    #[test]
    fn split_and_padding_test1() {
        let v = vec![1, 2, 3, 4];
        let r = split_and_padding(&v, 6);

        assert!(r.len() == 6);

        assert!(r[0].inner == &v[0..1]);
        assert!(r[1].inner == &v[1..2]);
        assert!(r[2].inner == &v[2..3]);
        assert!(r[3].inner == &v[3..4]);
        assert!(r[4].len == 0);
        assert!(r[5].len == 0);
    }

    #[test]
    fn split_and_padding_test2() {
        let v = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let r = split_and_padding(&v, 3);

        assert!(r.len() == 3);
        assert!(r[0].inner == &v[0..3]);
        assert!(r[1].inner == &v[3..6]);
        assert!(r[2].inner == &v[6..8]);
        assert!(r[2].len == 2);
    }
}
