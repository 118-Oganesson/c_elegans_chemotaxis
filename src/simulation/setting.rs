use rand::{rng, Rng};
use serde::{Deserialize, Serialize};

// シミュレーション設定を表す構造体Setting
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Setting {
    pub alpha: (f64, f64),    // 濃度勾配の範囲
    pub c_0: f64,             // ガウス濃度の設定
    pub lambda: f64,          // ガウス濃度の設定
    pub x_peak: f64,          // 勾配のピークのx座標 /cm
    pub y_peak: f64,          // 勾配のピークのy座標 /cm
    pub dt: f64,              // 時間刻みの幅 /s
    pub periodic_time: f64,   // 移動の1サイクルの継続時間 /s
    pub frequency: f64,       // 方向の平均周波数 /Hz
    pub velocity: f64,        // 線虫の速度 /cm/s
    pub simulation_time: f64, // シミュレーション時間 /s
    pub time_constant: f64,   // 時定数 /s
}

// 固定されたシミュレーション設定を表す構造体Const
pub struct Const {
    pub alpha: f64,           // 線形濃度の勾配
    pub c_0: f64,             // ガウス濃度の設定
    pub lambda: f64,          // ガウス濃度の設定
    pub x_peak: f64,          // 勾配のピークのx座標 /cm
    pub y_peak: f64,          // 勾配のピークのy座標 /cm
    pub dt: f64,              // 時間刻みの幅 /s
    pub periodic_time: f64,   // 移動の1サイクルの継続時間 /s
    pub frequency: f64,       // 方向の平均周波数 /Hz
    pub velocity: f64,        // 線虫の速度 /cm/s
    pub simulation_time: f64, // シミュレーション時間 /s
    pub time_constant: f64,   // 時定数 /s
}

impl Setting {
    // Setting構造体のデータを基にConst構造体を生成するメソッド
    pub fn const_new(&self) -> Const {
        // 乱数生成器を初期化
        let mut rng = rng();
        // alphaの範囲内で乱数を生成
        let a: f64 = rng.random_range(self.alpha.0..self.alpha.1);

        // Const構造体を生成して返す
        Const {
            alpha: a,
            c_0: self.c_0,
            lambda: self.lambda,
            x_peak: self.x_peak,
            y_peak: self.y_peak,
            dt: self.dt,
            periodic_time: self.periodic_time,
            frequency: self.frequency,
            velocity: self.velocity,
            simulation_time: self.simulation_time,
            time_constant: self.time_constant,
        }
    }
}
