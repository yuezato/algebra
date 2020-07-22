/// Finding Maximum Edges BiClique

use crate::bit_based_gf8::Command;
use std::collections::{HashMap, HashSet};
use rand::seq::SliceRandom;

#[derive(Debug)]
struct Data {
    // rows[to] = { from1, from2, ... }
    rows: HashMap<usize, HashSet<usize>>,

    // cols[from] = { to1, to2, ... }
    cols: HashMap<usize, HashSet<usize>>,
}

fn ops_len(data: &Data) -> usize {
    let mut acc = 0;

    for (_, froms) in &data.rows {
        acc += froms.len();
    }

    acc
}

fn ops_len2(data: &Data) -> usize {
    let mut acc = 0;

    for (_, tos) in &data.cols {
        acc += tos.len();
    }

    acc
}

fn update(ms: &mut HashMap<usize, HashSet<usize>>, key: usize, val: usize) {
    if let Some(set) = ms.get_mut(&key) {
        set.insert(val);
    } else {
        let mut set = HashSet::new();
        set.insert(val);
        ms.insert(key, set);
    }
}

fn commands_to_data(commands: &[Command]) -> Data {
    let mut rows = HashMap::new();
    let mut cols = HashMap::new();

    for c in commands {
        update(&mut rows, c.to(), c.from());
        update(&mut cols, c.from(), c.to());
    }
    
    Data { rows, cols }
}

fn data_to_commands(data: &Data) -> Vec<Command> {
    let mut commands = Vec::new();

    for (dst, srcs) in data.rows.iter() {
        let srcs: Vec<usize> = srcs.iter().cloned().collect();
        commands.push(Command::Copy(srcs[0], *dst));
        
        if srcs.len() > 1 {
            for src in &srcs[1..] {
                commands.push(Command::Xor(*src, *dst));
            }
        }
    }

    commands
}

fn expand_from_core(data: &Data, p: Vec<&HashSet<usize>>) -> (Vec<usize>, Vec<usize>) {
    let mut vecI: Vec<usize> = Vec::new();
    let mut vecJ: Vec<usize> = Vec::new();

    // println!("p = {:?}", p);
    
    // これintersectとった方が良くない?
    for col in data.cols.keys() {
        if p.iter().all(|set| set.contains(col)) {
            vecJ.push(*col);
        }
    }

    // println!("vecJ = {:?}", vecJ);

    // 非自明な場合だけ
    // この部分の処理をする
    if vecJ.len() > 1 {
        let vecJ_inner: Vec<&HashSet<usize>> =
            data.cols.iter()
            .filter(|(key, val)| vecJ.contains(key))
            .map(|(_, val)| val).collect();
        
        for row in data.rows.keys() {
            if vecJ_inner.iter().all(|set| set.contains(row)) {
                vecI.push(*row);
            }
        }
    }

    // println!("vecI = {:?}", vecI);
    
    (vecI, vecJ)
}

fn has_merit(nr_to: usize, nr_from: usize) -> bool {
    let naiive = nr_to * nr_from;
    
    // xored to buffer + broadcast to buffers
    // Notice: zeroed buffer is not needed.
    // in real, we do copy instead first xor
    let optimized = nr_from + nr_to;
    
    naiive > optimized
}

fn merit(nr_to: usize, nr_from: usize) -> usize {
    (nr_to * nr_from) - (nr_from + nr_to)
}

fn monte_carlo_maximum_edge_biclique(data: &Data) -> (Vec<usize>, Vec<usize>) {
    let mut ans_row = Vec::new();
    let mut ans_col = Vec::new();

    let rows: Vec<usize> = data.rows.keys().cloned().collect();
    let mut rng = rand::thread_rng();
    
    for _ in 0..30000 {
        for len in 2..=10 {
            let mut new_rows = Vec::new();

            for row in &rows {
                // 自明でないrowだけが選択対象になるようにする
                if data.rows[row].len() > 1 {
                    new_rows.push(*row);
                }
            }

            new_rows.shuffle(&mut rng);

            let rows = &new_rows[0..std::cmp::min(new_rows.len(), len)];
            
            let mut p: Vec<&HashSet<usize>> = Vec::new();
            for (key, value) in &data.rows {
                if rows.contains(key) {
                    p.push(value);
                }
            }

            // dump_data(&data);
            // println!("new_rows = {:?}", rows);
            
            let (row, col) = expand_from_core(data, p);
            if has_merit(row.len(), col.len()) &&
                merit(ans_row.len(), ans_col.len()) < merit(row.len(), col.len()) {
                ans_row = row;
                ans_col = col;
            }
        }
    }

    (ans_row, ans_col)
}

fn to_sorted_vec<'a>(s: impl Iterator<Item = &'a usize>) -> Vec<usize> {
    let mut v: Vec<usize> = s.cloned().collect();
    v.sort();
    v
}

fn dump_data(data: &Data) {
    let rows = to_sorted_vec(data.rows.keys());
            
    for to in &rows {
        if data.rows[to].len() > 1 {
            let value = to_sorted_vec(data.rows[to].iter());
            println!("{} <= {:?}", to, value);
        }
    }
}

fn generate_commands(dsts: &Vec<usize>, srcs: &Vec<usize>) -> Vec<Command> {
    let mut commands = Vec::new();
    
    // bufferをclearしてxoringするのではなく
    // 最初の一つは copy する
    commands.push(Command::Copy(srcs[0], 0xff));

    // 残りをxoringする
    for src in &srcs[1..] {
        commands.push(Command::Xor(*src, 0xff));
    }

    // bufferからbroadcastする
    for dst in dsts {
        commands.push(Command::Copy(0xff, *dst));
    }
    
    commands
}

fn shrink_data_by_MEB(data: &mut Data) -> Vec<Command> {
    use std::iter::FromIterator;

    let mut iter = 0;
    let mut true_delete = 0;
    let mut optimized_command = Vec::new();
    
    loop {
        println!("iter = {}", iter);
        println!("current_ops1 = {}, current_ops2 = {}", ops_len(&data), ops_len2(&data));
        // dump_data(&data);
        
        let (ans_row, ans_col) = monte_carlo_maximum_edge_biclique(&data);
        iter += 1;
        
        let row_set = HashSet::<usize>::from_iter(ans_row.iter().cloned());
        let col_set = HashSet::<usize>::from_iter(ans_col.iter().cloned());

        if has_merit(ans_row.len(), ans_col.len()) {
            true_delete += merit(ans_row.len(), ans_col.len());

            optimized_command.append(&mut generate_commands(&ans_row, &ans_col));
            
            println!("delete row.len() = {}, col.len() = {}", ans_row.len(), ans_col.len());
            
            for row in &ans_row {
                let set = data.rows.get_mut(row).unwrap();
                let subtracted = HashSet::<usize>::from_iter(set.difference(&col_set).cloned());
                *set = subtracted;
            }
            
            for col in &ans_col {
                let set = data.cols.get_mut(col).unwrap();
                let subtracted = HashSet::<usize>::from_iter(set.difference(&row_set).cloned());
                *set = subtracted;
            }
        } else {
            println!("total iter = {}, truly deleting = {}", iter, true_delete);
            dump_data(&data);
            return optimized_command;
        }
    }        
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::Iterator;
    use crate::vandermonde::*;
    use crate::fin_field::*;
    use crate::field::*;
    use crate::matrix::*;
    use crate::bit_based_gf8::*;

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
    fn maximum_edge_biclique_test1() {
        let m = make_vandermonde_matrix(10, 4);
        let commands = matrix_to_commands(&m);
        let data = commands_to_data(&commands);
        
        let (ans_row, ans_col) = monte_carlo_maximum_edge_biclique(&data);
        println!("score = {}, row = {}, col = {}",
                 ans_row.len() * ans_col.len(),
                 ans_row.len(), ans_col.len());
        
    }

    #[test]
    fn maximum_edge_biclique_test2() {
        let mut m = make_vandermonde_matrix(10, 4);
        m.drop_columns(vec![0, 1, 2, 3]);
        let inv = m.inverse().unwrap();
        
        let commands = matrix_to_commands(&inv);
        let data = commands_to_data(&commands);
        
        let (ans_row, ans_col) = monte_carlo_maximum_edge_biclique(&data);
        println!("commands.len() = {}", commands.len());
        println!("score = {}, row = {}, col = {}",
                 ans_row.len() * ans_col.len(),
                 ans_row.len(), ans_col.len());
        
    }

    #[test]
    fn shrink_data_by_meb_test1() {
        let mut m = make_vandermonde_matrix(10, 4);
        m.drop_columns(vec![0, 1, 2, 3]);
        let inv = m.inverse().unwrap();
        
        let commands = matrix_to_commands(&inv);
        let mut data = commands_to_data(&commands);

        println!("original #ops = {}, #ops = {}", commands.len(), ops_len(&data));
        let mut optimized = shrink_data_by_MEB(&mut data);
        println!("original #ops = {}, shrinked #ops = {}", commands.len(), ops_len(&data));
        optimized.append(&mut data_to_commands(&data));
        dbg!(optimized);
    }
    
    #[test]
    fn maximum_edge_biclique_test3() {
        let commands = vec![
            Command::Xor(0, 10),
            Command::Xor(1, 10),
            Command::Xor(2, 10),
            Command::Xor(0, 11),
            Command::Xor(1, 11),
            Command::Xor(2, 11),
            Command::Xor(0, 12),
            Command::Xor(1, 12),
            Command::Xor(2, 12),
            ];
        let data = commands_to_data(&commands);
        
        let (ans_row, ans_col) = monte_carlo_maximum_edge_biclique(&data);
        println!("commands.len() = {}", commands.len());
        println!("score = {}, row = {}, col = {}",
                 ans_row.len() * ans_col.len(),
                 ans_row.len(), ans_col.len());
    }   
    
    #[test]
    fn update_test1() {
        let mut ms = HashMap::new();

        update(&mut ms, 0, 1);

        let keys: Vec<usize> = to_sorted_vec(ms.keys());
        assert!(keys == vec![0]);

        let vals: Vec<usize> = to_sorted_vec(ms[&0].iter());
        assert!(vals == vec![1]);

        update(&mut ms, 0, 1);
        let keys: Vec<usize> = to_sorted_vec(ms.keys());
        assert!(keys == vec![0]);

        let vals: Vec<usize> = to_sorted_vec(ms[&0].iter());
        assert!(vals == vec![1]);
    }

    #[test]
    fn update_test2() {
        let mut ms = HashMap::new();

        update(&mut ms, 0, 1);

        let keys: Vec<usize> = to_sorted_vec(ms.keys());
        assert!(keys == vec![0]);

        let vals: Vec<usize> = to_sorted_vec(ms[&0].iter());
        assert!(vals == vec![1]);

        update(&mut ms, 0, 2);
        let keys: Vec<usize> = to_sorted_vec(ms.keys());
        assert!(keys == vec![0]);

        let vals: Vec<usize> = to_sorted_vec(ms[&0].iter());
        assert!(vals == vec![1, 2]);
    }

    #[test]
    fn update_test3() {
        let mut ms = HashMap::new();

        update(&mut ms, 0, 1);

        let keys: Vec<usize> = to_sorted_vec(ms.keys());
        assert!(keys == vec![0]);

        let vals: Vec<usize> = to_sorted_vec(ms[&0].iter());
        assert!(vals == vec![1]);

        update(&mut ms, 1, 3);
        let keys: Vec<usize> = to_sorted_vec(ms.keys());
        assert!(keys == vec![0, 1]);

        let vals0: Vec<usize> = to_sorted_vec(ms[&0].iter());
        assert!(vals0 == vec![1]);
        
        let vals1: Vec<usize> = to_sorted_vec(ms[&1].iter());
        assert!(vals1 == vec![3]);
    }

    #[test]
    fn commands_to_data_test1() {
        let mut commands = Vec::new();
        commands.push(Command::Xor(0, 1));
        commands.push(Command::Xor(2, 5));
        commands.push(Command::Xor(0, 2));
        commands.push(Command::Xor(3, 1));

        let data = commands_to_data(&commands);

        let tos = data.rows;
        let tos_keys: Vec<usize> = to_sorted_vec(tos.keys());
        assert!(tos_keys == vec![1, 2, 5]);

        assert!(to_sorted_vec(tos[&1].iter()) == vec![0, 3]);
        assert!(to_sorted_vec(tos[&2].iter()) == vec![0]);
        assert!(to_sorted_vec(tos[&5].iter()) == vec![2]);
        
        let froms = data.cols;
        let froms_keys: Vec<usize> = to_sorted_vec(froms.keys());
        assert!(froms_keys == vec![0, 2, 3]);
    }

    /*
     * 0 -> 1, 0 -> 2
     * 2 -> 5, 3 -> 1
     */
    fn data1() -> Data {
        let mut commands = Vec::new();
        commands.push(Command::Xor(0, 1));
        commands.push(Command::Xor(2, 5));
        commands.push(Command::Xor(0, 2));
        commands.push(Command::Xor(3, 1));

        commands_to_data(&commands)
    }
    
    #[test]
    fn expand_from_core_test1() {
        let data = data1();
        let p = vec![&data.rows[&1]];

        let (vI, vJ) = expand_from_core(&data, p);
        dbg!(vI);
        dbg!(vJ);
    }
}
