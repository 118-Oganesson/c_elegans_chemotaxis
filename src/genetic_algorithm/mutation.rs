use crate::genetic_algorithm::ga::*;
use crate::simulation::*;
use rand::Rng;
use rand_distr::{Distribution, Normal};

// mutation関数は遺伝子に突然変異を起こします
pub fn mutation(gene: &Ga, rate: f64, mean_std: (f64, f64)) -> Ga {
    let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
    let normal = Normal::new(mean_std.0, mean_std.1).expect("Failed to create normal distribution");

    let mut mutated_gene: Gene = gene.gene.clone();

    for val in mutated_gene.gene.iter_mut() {
        if rng.gen::<f64>() < rate {
            let delta: f64 = normal.sample(&mut rng);
            *val += delta;
            *val = (*val).clamp(-1.0, 1.0);
        }
    }

    Ga {
        value: 0.0,
        gene: mutated_gene,
    }
}

// mutation_aiy_aiz_negative関数はAIY-AIZシナプスの遺伝子を抑制に固定した上で突然変異を起こします
pub fn mutation_aiy_aiz_negative(gene: &Ga, rate: f64, mean_std: (f64, f64)) -> Ga {
    let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
    let normal: Normal<f64> =
        Normal::new(mean_std.0, mean_std.1).expect("Failed to create normal distribution");

    let mut mutated_gene: Gene = gene.gene.clone();

    for (index, val) in mutated_gene.gene.iter_mut().enumerate() {
        if index > 11 && index < 14 {
            // 遺伝子番号12〜13(w_02, w_13の結合)
            if rng.gen::<f64>() < rate {
                let delta: f64 = normal.sample(&mut rng);
                *val += delta;
                if *val > 0.0 {
                    *val -= delta
                } else if *val < -1.0 {
                    *val = -1.0
                }
            }
        } else if rng.gen::<f64>() < rate {
            let delta: f64 = normal.sample(&mut rng);
            *val += delta;
            *val = (*val).clamp(-1.0, 1.0);
        }
    }

    Ga {
        value: 0.0,
        gene: mutated_gene,
    }
}
