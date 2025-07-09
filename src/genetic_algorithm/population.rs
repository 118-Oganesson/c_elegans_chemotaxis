use crate::genetic_algorithm::ga::*;
use crate::simulation::*;
use rand::{rng, Rng};

// population_new関数はランダムな遺伝子を持つ個体群を生成します
pub fn population_new(gen_size: usize, pop_size: usize) -> Vec<Ga> {
    let mut rng = rng();
    let population: Vec<Ga> = (0..pop_size)
        .map(|_| {
            let gene: Gene = (0..gen_size)
                .map(|_| rng.random_range(-1.0..1.0))
                .collect::<Gene>();
            Ga { value: 0.0, gene }
        })
        .collect();
    population
}

// population_new_aiy_aiz_negative関数はAIY-AIZシナプスの遺伝子を抑制に固定する個体群を生成します
pub fn population_new_aiy_aiz_negative(gen_size: usize, pop_size: usize) -> Vec<Ga> {
    let mut rng = rng();
    let population: Vec<Ga> = (0..pop_size)
        .map(|_| {
            let mut gene: Gene = (0..gen_size)
                .map(|_| rng.random_range(-1.0..1.0))
                .collect::<Gene>();

            // 遺伝子番号12〜13(w_02, w_13の結合)をマイナスに設定
            for i in 12..14 {
                if gene.gene[i] > 0.0 {
                    gene.gene[i] = -gene.gene[i];
                }
            }

            Ga { value: 0.0, gene }
        })
        .collect();
    population
}
