use crate::simulation::*;
use serde::{Deserialize, Serialize};

// Ga構造体は遺伝的アルゴリズムの中での個体を表します
#[derive(Clone, Debug)]
pub struct Ga {
    pub value: f64, // 個体の評価値(CI値)
    pub gene: Gene, // 遺伝子
}

// Gajson構造体は遺伝的アルゴリズムの中での個体をJSON形式で表します
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gajson {
    pub value: f64,     // 個体の評価値(CI値)
    pub gene: Vec<f64>, // 遺伝子を表す配列
}
