use crate::genetic_algorithm::ga::*;
use crate::simulation::*;
use rand::{rng, Rng};

// two_point_crossover関数は2点交叉を実行し、2つの子個体を生成します
pub fn two_point_crossover(parent1: &Ga, parent2: &Ga) -> Vec<Ga> {
    // ランダムな交叉点を決定
    let mut rng = rng();
    let crossover_points: (usize, usize) = (
        rng.random_range(0..parent1.gene.gene.len()),
        rng.random_range(0..parent1.gene.gene.len()),
    );
    let (start, end) = if crossover_points.0 < crossover_points.1 {
        (crossover_points.0, crossover_points.1)
    } else {
        (crossover_points.1, crossover_points.0)
    };

    // 2つの子個体の遺伝子を生成
    let mut child_gene1: Gene = parent1.gene.clone();
    let mut child_gene2: Gene = parent2.gene.clone();
    for i in start..end {
        child_gene1.gene[i] = parent2.gene.gene[i];
        child_gene2.gene[i] = parent1.gene.gene[i];
    }

    // 子個体をGa構造体に格納して返す
    let child1: Ga = Ga {
        value: 0.0,
        gene: child_gene1,
    };
    let child2: Ga = Ga {
        value: 0.0,
        gene: child_gene2,
    };
    vec![child1, child2]
}
