use crate::analysis::bearing::*;
use crate::analysis::curving_rate::*;
use crate::analysis::normal_gradient::*;
use crate::analysis::translational_gradient::*;
use crate::simulation::*;
use rand::{rng, Rng};
use std::collections::VecDeque;
use std::f64::consts::PI;

#[allow(clippy::too_many_arguments)]
pub fn klinotaxis_analysis(
    gene: &Gene,
    setting: &Setting,
    mode: usize,
    periodic_number: usize,
    periodic_number_drain: usize,
    delta: f64,
    analysis_use_function: &i32,
) -> Vec<Vec<f64>> {
    //定数
    let constant: Const = setting.const_new();
    //遺伝子の受け渡し
    let weight: GeneConst = gene.scaling();
    //時間に関する定数をステップ数に変換
    let time: Time = time_new(&weight, &constant);

    //配列の宣言
    let mut y: [[f64; 8]; 2] = [[0.0; 8]; 2];
    let mut mu: [f64; 2] = [0.0; 2];
    let mut phi: [f64; 2] = [0.0; 2];
    let mut r: [[f64; 2]; 2] = [[0.0; 2]; 2];

    //Vecの宣言
    let mut r_vec: Vec<[f64; 2]> = vec![[0.0; 2]; 1];
    let mut mu_vec: Vec<f64> = vec![0.0];

    // concentration関数の選択
    let concentration_fn: fn(&Const, f64, f64) -> f64 = match mode {
        0 => concentration,
        1 => gauss_concentration,
        2 => two_gauss_concentration,
        _ => panic!("Invalid mode: {mode}"),
    };

    //初期濃度の履歴生成
    let mut c_vec: VecDeque<f64> = VecDeque::new();
    for _ in 0..time.n_time + time.m_time {
        c_vec.push_back(concentration_fn(&constant, 0.0, 0.0));
    }

    //運動ニューロンの初期活性を0～1の範囲でランダム化
    let mut rng_init = rng();
    rng_init.fill(&mut y[0][4..]);

    //ランダムな向きで配置
    let mut rng_dir = rng();
    mu[0] = rng_dir.random_range(0.0..2.0 * PI);

    //オイラー法
    for k in 0..time.simulation_time - 1 {
        //濃度の更新
        c_vec.pop_front();
        c_vec.push_back(concentration_fn(&constant, r[0][0], r[0][1]));

        let y_on_off: [f64; 2] = y_on_off(&constant, &weight, &time, &c_vec);
        let y_osc: f64 = y_osc(&constant, k as f64 * constant.dt);

        for i in 0..8 {
            let mut synapse: f64 = 0.0;
            let mut gap: f64 = 0.0;
            for j in 0..8 {
                synapse += weight.w[j][i] * sigmoid(y[0][j] + weight.theta[j]);
                gap += weight.g[j][i] * (y[0][j] - y[0][i]);
            }
            //外部からの入力
            let input: f64 = weight.w_on[i] * y_on_off[0]
                + weight.w_off[i] * y_on_off[1]
                + weight.w_osc[i] * y_osc;
            //ニューロンの膜電位の更新
            y[1][i] =
                y[0][i] + (-y[0][i] + synapse + gap + input) / constant.time_constant * constant.dt;
        }

        //方向の更新
        let d: f64 = sigmoid(y[0][5] + weight.theta[5]) + sigmoid(y[0][6] + weight.theta[6]);
        let v: f64 = sigmoid(y[0][4] + weight.theta[4]) + sigmoid(y[0][7] + weight.theta[7]);
        phi[1] = phi[0];
        phi[0] = weight.w_nmj * (d - v);
        mu[1] = mu[0] + phi[0] * constant.dt;

        //位置の更新
        r[1][0] = r[0][0] + constant.velocity * (mu[0]).cos() * constant.dt;
        r[1][1] = r[0][1] + constant.velocity * (mu[0]).sin() * constant.dt;

        //Vecへの追加
        r_vec.push(r[1]);
        mu_vec.push(mu[1]);

        //更新
        for i in 0..8 {
            y[0][i] = y[1][i];
        }
        mu[0] = mu[1];
        for i in 0..2 {
            r[0][i] = r[1][i];
        }
    }

    let mut x_axis_value: Vec<f64> = Vec::new();
    if *analysis_use_function == 0 {
        // Bearing
        x_axis_value = bearing(&r_vec, &constant, &time, &periodic_number);
    } else if *analysis_use_function == 1 {
        // Normal gradient
        x_axis_value = normal_gradient(&r_vec, &constant, &time, &periodic_number, delta, mode);
    } else if *analysis_use_function == 2 {
        // Translational gradient
        x_axis_value =
            translational_gradient(&r_vec, &constant, &time, &periodic_number, delta, mode);
    }

    // Turning bias
    let mut y_axis_value: Vec<f64> = curving_rate_bear(&r_vec, &time, &periodic_number);

    // 先頭の要素を削除
    x_axis_value.drain(..periodic_number_drain * time.periodic_time);
    y_axis_value.drain(..periodic_number_drain * time.periodic_time);

    // 結果
    let result: Vec<Vec<f64>> = vec![x_axis_value, y_axis_value];

    result
}
