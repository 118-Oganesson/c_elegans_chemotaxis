use serde::{Deserialize, Serialize};

// Gasetting構造体は遺伝的アルゴリズムの設定を表します
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gasetting {
    pub average: i32,    // CI測定時の平均回数
    pub gen_size: usize, // 遺伝子の長さ
    pub ga_count: usize, // 遺伝的アルゴリズムの実行回数
    pub n_gen: usize,    // 世代数
    pub pop_size: usize, // 個体数
    pub sel_top: usize,  // 選択する上位個体の数
    pub mat_pb: f64,     // 交差する確率
    pub mut_pb: f64,     // 変異する確率
    pub re_val: usize,   // 全個体を再評価する間隔
}
