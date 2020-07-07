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

    println!("before: \n{}", generator.matrix().dump());

    let m = decode_matrix(generator.take_matrix(), &alive).unwrap();

    println!("after: \n{}", m.dump());

    let mul_table = &table[&alive];

    mom2(&m, mul_table, &encoded_data).into_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

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
