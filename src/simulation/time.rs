use crate::simulation::gene::*;
use crate::simulation::setting::*;

// 構造体Timeはシミュレーションの時間関連の設定を保持するための構造体
#[allow(dead_code)]
pub struct Time {
    pub n_time: usize, // 感覚ニューロンの時間スケールに対応するタイムステップ数
    pub m_time: usize, // 運動ニューロンの時間スケールに対応するタイムステップ数
    pub simulation_time: usize, // シミュレーション全体のタイムステップ数
    pub f_inv_time: usize, // 周波数の逆数に対応するタイムステップ数
    pub periodic_time: usize, // 移動の1サイクルのタイムステップ数
}

// 関数time_newはGeneConstとConstの設定を基にTime構造体を生成する
pub fn time_new(weight: &GeneConst, constant: &Const) -> Time {
    let n_time: usize = (weight.n / constant.dt).floor() as usize;
    let m_time: usize = (weight.m / constant.dt).floor() as usize;
    let simulation_time: usize = (constant.simulation_time / constant.dt).floor() as usize;
    let f_inv_time: usize = (1.0 / constant.frequency / constant.dt).floor() as usize;
    let periodic_time: usize = (constant.periodic_time / constant.dt).floor() as usize;

    Time {
        n_time,
        m_time,
        simulation_time,
        f_inv_time,
        periodic_time,
    }
}
