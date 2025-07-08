use crate::genetic_algorithm::ga::*;
use crate::simulation::*;
use rayon::prelude::*;

impl Ga {
    // fitnessメソッドは個体の適応度を計算します
    pub fn fitness(&self, setting: &Setting, average: i32, version: i32) -> f64 {
        let mut sum: f64 = 0.0;
        // バージョン0の場合: 通常の化学性走化性指数を計算
        if version == 0 {
            for _ in 0..average {
                sum += chemotaxis_index(&self.gene, setting);
            }
        // バージョン1の場合: 波形チェックを含む化学性走化性指数を計算
        } else if version == 1 {
            for _ in 0..average {
                sum += chemotaxis_index_wave_check(&self.gene, setting);
            }
        }
        sum / average as f64
    }
}

// evaluate_fitness関数は個体群の適応度を評価します
pub fn evaluate_fitness(
    population: &mut Vec<Ga>,
    setting: &Setting,
    average: i32,
    version: i32,
) -> Vec<Ga> {
    population
        .par_iter_mut() // パラレルイテレータを使用して並列処理を行う
        .map(|ind| Ga {
            value: ind.fitness(setting, average, version),
            gene: ind.gene.clone(),
        })
        .collect()
}
