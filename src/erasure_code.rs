use crate::fin_field::*;
use crate::matrix::*;
use crate::reed_solomon::*;
use std::collections::HashMap;
use std::string::ToString;

/*
 * 1. alive情報から復号行列を作る
 * 2. alive情報から復号乗法表を作る
 */

#[derive(Clone, Hash, PartialEq, Eq)]
pub struct AliveBlocks {
    inner: u128,
    blocks: usize,
}

impl AliveBlocks {
    pub fn to_boolvec(&self) -> Vec<bool> {
        let mut v = vec![false; self.blocks];
        for (i, e) in v.iter_mut().enumerate() {
            *e = self.at(i);
        }
        v
    }

    pub fn from_boolvec(v: &[bool]) -> Self {
        let mut inner = 0;
        for (i, e) in v.iter().enumerate() {
            if *e {
                inner |= 1 << i;
            }
        }
        Self {
            inner,
            blocks: v.len(),
        }
    }

    pub fn from_alive_vec(blocks: usize, v: &[usize]) -> Self {
        let mut inner = 0;
        for e in v {
            inner |= 1 << *e;
        }
        Self { inner, blocks }
    }

    pub fn new(blocks: usize) -> Self {
        Self { inner: 0, blocks }
    }

    pub fn next(&self) -> Option<AliveBlocks> {
        let v = self.inner + 1;
        if v >= (1 << self.blocks) {
            None
        } else {
            Some(Self {
                inner: v,
                blocks: self.blocks,
            })
        }
    }

    pub fn blocks(&self) -> usize {
        self.blocks
    }

    pub fn num_livings(&self) -> u32 {
        self.inner.count_ones()
    }

    pub fn at(&self, idx: usize) -> bool {
        ((self.inner >> idx) & 1) == 1
    }

    pub fn at_mut(&mut self, idx: usize, b: bool) {
        if b {
            self.inner |= 1 << idx;
        } else {
            self.inner |= 1 << idx;
            self.inner ^= 1 << idx;
        }
    }

    pub fn alive_blocks(&self) -> Vec<usize> {
        let mut v = Vec::new();
        for i in 0..self.blocks {
            if self.at(i) {
                v.push(i);
            }
        }
        v
    }

    pub fn erased_blocks(&self) -> Vec<usize> {
        let mut v = Vec::new();
        for i in 0..self.blocks {
            if !self.at(i) {
                v.push(i);
            }
        }
        v
    }
}

pub fn decode_matrix<F: FiniteField + ToString>(
    mut generator: Matrix<F>,
    alive: &AliveBlocks,
) -> Option<Matrix<F>> {
    debug_assert!(alive.blocks() == generator.height());

    // width以上残っていなければ復元できない
    if (alive.num_livings() as usize) < generator.width() {
        return None;
    }

    // 前方優先にするのはできるだけ単位行列に近づけたいから
    let mut remove_blocks: Vec<usize> = alive.alive_blocks()[generator.width()..].to_vec();
    remove_blocks.append(&mut alive.erased_blocks());
    generator.drop_columns(remove_blocks);

    generator.inverse()
}

pub fn make_decode_multable<F: FiniteField + ToString>(
    generator: &Matrix<F>,
) -> HashMap<AliveBlocks, MulTable<F>> {
    let blocks = generator.height();
    let mut alive = AliveBlocks::new(blocks);

    let mut memo = HashMap::new();

    loop {
        if alive.num_livings() as usize >= generator.width() {
            if let Some(m) = decode_matrix(generator.clone(), &alive) {
                memo.insert(alive.clone(), MulTable::build(&m));
            }
        }
        if let Some(alive_) = alive.next() {
            alive = alive_;
        } else {
            break memo;
        }
    }
}

fn dot_prod_row_and_matrix<F: FiniteField>(v: &[F], m: &[&[u8]]) -> Vec<u8> {
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

#[allow(non_snake_case)]
pub fn decode_by_RSV2<F: FiniteField + ToString>(
    generator: Generator<F>,
    data: Vec<Encoded>,
) -> Vec<u8> {
    let alive_nums: Vec<usize> = data.iter().map(|e| e.block_num()).collect();
    let block_nums = generator.matrix().height();
    let alive = AliveBlocks::from_alive_vec(block_nums, &alive_nums);

    let mut encoded_data: Vec<&[u8]> = data.iter().map(|e| e.data()).collect();

    let generator: Matrix<F> = generator.take_matrix();

    let erased_data_blocks: Vec<usize> = (0..generator.width()).filter(|x| !alive.at(*x)).collect();

    /*
     * データブロックの消失が0なら元データが事実上残っているので
     * それを返してやる
     */
    if erased_data_blocks.is_empty() {
        return encoded_data[0..generator.width()].concat();
    }

    /*
     * データブロックの消失が1未満 かつ
     * parity先頭が消えてない場合は最適化が可能
     */
    // データブロックの消失数
    if erased_data_blocks.len() == 1 && alive.at(generator.width()) {
        let reconstruct_block = erased_data_blocks[0];

        // 全部1のハズなのでassertを書くこと
        let parity_top_row: &Vec<F> = generator.column_vec(generator.width()).as_vec();

        let data: &[&[u8]] = &encoded_data[0..generator.width()];
        let reconstructed: Vec<u8> = dot_prod_row_and_matrix(parity_top_row, data);

        debug_assert!(reconstructed.len() == encoded_data[0].len());

        encoded_data.insert(reconstruct_block, &reconstructed);
        return encoded_data[0..generator.width()].concat();
    }

    /*
     * それ以外の場合は逆行列を求めて計算する
     */
    let m = decode_matrix(generator, &alive).unwrap();
    mom(&m, &encoded_data).into_vec()
}

#[allow(non_snake_case)]
pub fn decode_by_RSV3<F: FiniteField + ToString>(
    generator: Generator<F>,
    data: Vec<Encoded>,
) -> Vec<u8> {
    let generator: Matrix<F> = generator.take_matrix();
    let nr_data = generator.width();
    let block_nums = generator.height();
    let parity_top_row: Vec<F> = generator.column_vec(generator.width()).as_vec().clone();
    
    let alive_nums: Vec<usize> = data.iter().map(|e| e.block_num()).collect();
    let alive = AliveBlocks::from_alive_vec(block_nums, &alive_nums);

    let encoded_data: Vec<&[u8]> = data.iter().map(|e| e.data()).collect();

    // FIX: Vec<Vec<u8>>は回避したい
    let mut result_data: Vec<Vec<u8>> = Vec::new();
    for d in &data {
        if d.block_num() < nr_data {
            result_data.push(d.data().to_vec());
        }
    }
    
    let mut erased_data_blocks: Vec<usize> = (0..generator.width()).filter(|x| !alive.at(*x)).collect();

    /*
     * データブロックの消失が0なら元データが事実上残っているので
     * それを返してやる
     */
    if erased_data_blocks.is_empty() {
        return encoded_data[0..generator.width()].concat();
    }

    let decoder = decode_matrix(generator, &alive).unwrap();

    // parity先頭が存在している場合
    if alive.at(nr_data) {
        let topmost_parity_block: &[u8] = data.iter().find(|x| x.block_num() == nr_data).unwrap().data();
        result_data.push(topmost_parity_block.to_vec());
        for i in 0..nr_data-1 {
            if !alive.at(i) {
                // 消失しているのでデコードして適切な位置に入れる
                let row = decoder.column_vec(i).as_vec();

                let decoded_row: Vec<u8>= dot_prod_row_and_matrix(row, &encoded_data);
                result_data.insert(i, decoded_row);

                let _ = erased_data_blocks.remove(0);
            }
            if erased_data_blocks.len() == 1 {
                // 残り一つになったので topmost parity を使った復号を行う
                // この段階では他の全ての元データが揃っていて欲しい
                let reconstruct_block = erased_data_blocks[0];
                let ref_data: Vec<&[u8]> = result_data.iter().map(|v| &v[..]).collect();
                let reconstructed: Vec<u8> = dot_prod_row_and_matrix(&parity_top_row, &ref_data);
                debug_assert!(reconstructed.len() == result_data[0].len());
                result_data.insert(reconstruct_block, reconstructed);
                return result_data[0..nr_data].concat();
            }
        }
        unreachable!("unreachable");
    } else {
        mom(&decoder, &encoded_data).into_vec()
    }
}

/*
 * generator: 生成行列、符号化に使った行列そのまま
 */
pub fn decode_by_table<F: FiniteField + ToString>(
    generator: Generator<F>,
    table: &HashMap<AliveBlocks, MulTable<F>>,
    data: Vec<Encoded>,
) -> Vec<u8> {
    let alive_nums: Vec<usize> = data.iter().map(|e| e.block_num()).collect();
    let block_nums = generator.matrix().height();
    let alive = AliveBlocks::from_alive_vec(block_nums, &alive_nums);

    let encoded_data: Vec<&[u8]> = data.iter().map(|e| e.data()).collect();

    // println!("before: \n{}", generator.matrix().dump());

    let m = decode_matrix(generator.take_matrix(), &alive).unwrap();

    // println!("after: \n{}", m.dump());

    let mul_table = &table[&alive];

    mom2(&m, mul_table, &encoded_data).into_vec()
}

#[allow(non_snake_case)]
pub fn decode_by_table2<F: FiniteField + ToString>(
    generator: Generator<F>,
    table: &HashMap<AliveBlocks, MulTable<F>>,
    data: Vec<Encoded>,
) -> Vec<u8> {
    let alive_nums: Vec<usize> = data.iter().map(|e| e.block_num()).collect();
    let block_nums = generator.matrix().height();
    let alive = AliveBlocks::from_alive_vec(block_nums, &alive_nums);

    let mut encoded_data: Vec<&[u8]> = data.iter().map(|e| e.data()).collect();

    let generator: Matrix<F> = generator.take_matrix();

    let erased_data_blocks: Vec<usize> = (0..generator.width()).filter(|x| !alive.at(*x)).collect();

    if erased_data_blocks.is_empty() {
        return encoded_data[0..generator.width()].concat();
    }

    if erased_data_blocks.len() == 1 && alive.at(generator.width()) {
        let reconstruct_block = erased_data_blocks[0];

        // 全部1のハズなのでassertを書くこと
        let parity_top_row: &Vec<F> = generator.column_vec(generator.width()).as_vec();

        let data: &[&[u8]] = &encoded_data[0..generator.width()];
        let reconstructed: Vec<u8> = dot_prod_row_and_matrix(parity_top_row, data);

        debug_assert!(reconstructed.len() == encoded_data[0].len());

        encoded_data.insert(reconstruct_block, &reconstructed);
        return encoded_data[0..generator.width()].concat();
    }

    let m = decode_matrix(generator, &alive).unwrap();
    let mul_table = &table[&alive];

    mom2(&m, mul_table, &encoded_data).into_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1_decode_by_rsv2() {
        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();

        fn testfunc<F: FiniteField + HasPrimitiveElement + ToString>() {
            let data_size = 4;
            let parity_size = 2;

            const MB: usize = 1024 * 1024;
            let original_data: Vec<u8> = (0..1 * MB).map(|_| rand::random::<u8>()).collect();
            let (generator, encoded) = encode_by_RSV::<F>(data_size, parity_size, &original_data);

            // 何も消えていない状態での復号化
            let decoded: Vec<u8> = decode_by_RSV2::<F>(generator, encoded);
            assert_eq!(original_data, decoded);
        }
    }

    #[test]
    fn test2_decode_by_rsv2() {
        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();

        fn testfunc<F: FiniteField + HasPrimitiveElement + ToString>() {
            let data_size = 4;
            let parity_size = 2;

            const MB: usize = 1024 * 1024;
            let original_data: Vec<u8> = (0..1 * MB).map(|_| rand::random::<u8>()).collect();
            let (generator, mut encoded) =
                encode_by_RSV::<F>(data_size, parity_size, &original_data);

            // パリティブロック2つ目と
            // データブロック先頭を消す
            encoded.remove(5);
            encoded.remove(1);

            // データブロックが1つ消えているが
            // パリティブロック先頭は残っている場合のテスト
            let decoded: Vec<u8> = decode_by_RSV2::<F>(generator, encoded);
            assert!(original_data == decoded);
        }
    }

    #[test]
    fn test3_decode_by_rsv2() {
        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();

        fn testfunc<F: FiniteField + HasPrimitiveElement + ToString>() {
            let data_size = 4;
            let parity_size = 2;

            const MB: usize = 1024 * 1024;
            let original_data: Vec<u8> = (0..1 * MB).map(|_| rand::random::<u8>()).collect();
            let (generator, mut encoded) =
                encode_by_RSV::<F>(data_size, parity_size, &original_data);

            // パリティブロック先頭と
            // データブロック先頭を消す
            encoded.remove(5);
            encoded.remove(1);

            // 逆行列を使って計算する場合
            let decoded: Vec<u8> = decode_by_RSV2::<F>(generator, encoded);
            assert!(original_data == decoded);
        }
    }

        #[test]
    fn test1_decode_by_rsv3() {
        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();

        fn testfunc<F: FiniteField + HasPrimitiveElement + ToString>() {
            let data_size = 4;
            let parity_size = 2;

            const MB: usize = 1024 * 1024;
            let original_data: Vec<u8> = (0..1 * MB).map(|_| rand::random::<u8>()).collect();
            let (generator, encoded) = encode_by_RSV::<F>(data_size, parity_size, &original_data);

            // 何も消えていない状態での復号化
            let decoded: Vec<u8> = decode_by_RSV3::<F>(generator, encoded);
            assert_eq!(original_data, decoded);
        }
    }

    #[test]
    fn test2_decode_by_rsv3() {
        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();

        fn testfunc<F: FiniteField + HasPrimitiveElement + ToString>() {
            let data_size = 4;
            let parity_size = 2;

            const MB: usize = 1024 * 1024;
            let original_data: Vec<u8> = (0..1 * MB).map(|_| rand::random::<u8>()).collect();
            let (generator, mut encoded) =
                encode_by_RSV::<F>(data_size, parity_size, &original_data);

            // パリティブロック2つ目と
            // データブロック先頭を消す
            encoded.remove(5);
            encoded.remove(1);

            // データブロックが1つ消えているが
            // パリティブロック先頭は残っている場合のテスト
            let decoded: Vec<u8> = decode_by_RSV3::<F>(generator, encoded);
            assert!(original_data == decoded);
        }
    }

    #[test]
    fn test3_decode_by_rsv3() {
        testfunc::<GF_2_8>();
        testfunc::<GF_2_16_Val>();

        fn testfunc<F: FiniteField + HasPrimitiveElement + ToString>() {
            let data_size = 4;
            let parity_size = 2;

            const MB: usize = 1024 * 1024;
            let original_data: Vec<u8> = (0..1 * MB).map(|_| rand::random::<u8>()).collect();
            let (generator, mut encoded) =
                encode_by_RSV::<F>(data_size, parity_size, &original_data);

            // パリティブロック先頭と
            // データブロック先頭を消す
            encoded.remove(5);
            encoded.remove(1);

            // 逆行列を使って計算する場合
            let decoded: Vec<u8> = decode_by_RSV3::<F>(generator, encoded);
            assert!(original_data == decoded);
        }
    }

    #[test]
    fn alive_blocks_test1() {
        let f = false;
        let t = true;
        let mut ab = AliveBlocks::new(3);
        assert_eq!(ab.to_boolvec(), vec![f, f, f]);
        assert_eq!(ab.num_livings(), 0);

        ab = ab.next().unwrap();
        assert_eq!(ab.to_boolvec(), vec![t, f, f]);
        assert_eq!(ab.num_livings(), 1);

        ab = ab.next().unwrap();
        assert_eq!(ab.to_boolvec(), vec![f, t, f]);
        assert_eq!(ab.num_livings(), 1);

        ab = ab.next().unwrap();
        assert_eq!(ab.to_boolvec(), vec![t, t, f]);
        assert_eq!(ab.num_livings(), 2);

        ab = ab.next().unwrap();
        assert_eq!(ab.to_boolvec(), vec![f, f, t]);
        assert_eq!(ab.num_livings(), 1);

        ab = ab.next().unwrap();
        assert_eq!(ab.to_boolvec(), vec![t, f, t]);
        assert_eq!(ab.num_livings(), 2);

        ab = ab.next().unwrap();
        assert_eq!(ab.to_boolvec(), vec![f, t, t]);
        assert_eq!(ab.num_livings(), 2);

        ab = ab.next().unwrap();
        assert_eq!(ab.to_boolvec(), vec![t, t, t]);
        assert_eq!(ab.num_livings(), 3);

        assert!(ab.next().is_none());
    }

    #[test]
    fn alive_blocks_test2() {
        let f = false;
        let t = true;
        let mut ab = AliveBlocks::new(3);
        assert_eq!(ab.to_boolvec(), vec![f, f, f]);

        ab.at_mut(0, true);
        assert_eq!(ab.to_boolvec(), vec![t, f, f]);
        ab.at_mut(0, false);
        assert_eq!(ab.to_boolvec(), vec![f, f, f]);

        let mut ab = AliveBlocks::from_boolvec(&vec![f, t, t]);
        assert_eq!(ab.to_boolvec(), vec![f, t, t]);

        ab.at_mut(2, false);
        assert_eq!(ab.to_boolvec(), vec![f, t, f]);
    }
}
