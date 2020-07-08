use crate::univariate_polynomial::*;
use crate::field::*;
use crate::fin_field::*;

#[allow(non_camel_case_types)]
pub struct Bit_GF_2_8_impl {
    ppoly: Poly<GF_2>,
    table: Vec<[u8; 8]>,
}

fn xors(x: u8) -> u8 {
    (x.count_ones() % 2) as u8
}

impl Bit_GF_2_8_impl {
    pub fn new(ppoly: Poly<GF_2>) -> Self {
        Self {
            ppoly,
            table: Vec::new(),
        }
    }

    pub fn setup(&mut self) {
        for i in 0u8..=255u8 {
            let t = self.conv(i.into());
            self.table.push(t);
        }
    }
    
    pub fn mul(&self, p: GF_2_8, q: GF_2_8) -> GF_2_8 {
        // p[deg] = c_of_deg7 c_of_deg6 ... c_of_deg1 c_of_deg0
        let p: [u8; 8] = self.table[u8::from(p) as usize];

        // (q >> deg) = 係数
        let q: u8 = q.into();
        
        let mut r: u8 = 0;
        
        for deg in 0..8 {
            let v: u8 = p[deg] & q;
            let v: u8 = xors(v);
            // 以上2つで内積計算をしている
            // v == 0 or v == 1

            r |= v << deg;
        }

        r.into()
    }

    // matrix積の布石
    pub fn mul2(&self, p: GF_2_8, q: GF_2_8) -> GF_2_8 {
        let p: [u8; 8] = self.table[u8::from(p) as usize];
        let mut r: u8 = 0;

        for i in 0..8 {
            let row = p[i];

            for j in 0..8 {
                if (row >> j) & 1 == 1 {
                    r ^= ((u8::from(q) >> j) & 1) << i;
                } else {
                    // 0かけてxorなので何もしなくて良い
                }
            }
        }

        r.into()
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
        let mut r: [u8; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

        // 次数iの計算
        // assume: rv = &r[i]
        for (i, rv) in r.iter_mut().enumerate() {
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
            *rv = v;
        }

        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;
    
    fn ppoly() -> Poly<GF_2> {
        Poly::from_vec(vec![
            (8, GF_2::ONE), (4, GF_2::ONE), (3, GF_2::ONE), (2, GF_2::ONE), (0, GF_2::ONE)
        ])
    }
    
    #[test]
    fn mul_test() {
        let i1 = GF_2_8_impl::new(ppoly());
        let mut i2 = Bit_GF_2_8_impl::new(ppoly());
        i2.setup();

        for i in 0u8..=255u8 {
            for j in 0u8..=255u8 {
                let r1 = i1.mul(i.into(), j.into());
                let r2 = i2.mul(i.into(), j.into());
                let r3 = i2.mul2(i.into(), j.into());
                
                assert!(r1 == r2, "{} * {}", i, j);
                assert!(r2 == r3, "{} * {}", i, j);
            }
        }
    }

    #[test]
    fn speed_test() {
        let i1 = GF_2_8_impl::new(ppoly());
        let mut i2 = Bit_GF_2_8_impl::new(ppoly());
        i2.setup();

        let mut r1: u8 = 0;
        let mut r2: u8 = 0;
        let mut r3: u8 = 0;
        
        let t1 = Instant::now();

        for _ in 0..1000 {
            for i in 0u8..=255u8 {
                for j in 0u8..=255u8 {
                    r1 ^= u8::from(i1.mul(i.into(), j.into()));
                }
            }
        }
        println!("naiive mul = {:?}", t1.elapsed());

        let t2 = Instant::now();
        for _ in 0..1000 {
            for i in 0u8..=255u8 {
                for j in 0u8..=255u8 {
                    r2 ^= u8::from(i2.mul(i.into(), j.into()));
                }
            }
        }
        println!("bitmatrix mul = {:?}", t2.elapsed());

        let t3 = Instant::now();
        for _ in 0..1000 {
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
    
