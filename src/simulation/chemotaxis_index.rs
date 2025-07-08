use crate::simulation::function::*;
use crate::simulation::gene::*;
use crate::simulation::setting::*;
use crate::simulation::time::*;
use rand::{thread_rng, Rng};
use std::collections::VecDeque;
use std::f64::consts::PI;

// chemotaxis_index関数は化学的走化性指数を計算します
pub fn chemotaxis_index(gene: &Gene, setting: &Setting) -> f64 {
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

    //初期濃度の履歴生成
    let mut c_vec: VecDeque<f64> = VecDeque::new();
    for _ in 0..time.n_time + time.m_time {
        c_vec.push_back(concentration(&constant, 0.0, 0.0));
    }
    //運動ニューロンの初期活性を0～1の範囲でランダム化
    let _ = thread_rng().try_fill(&mut y[0][4..]);

    //ランダムな向きで配置
    let mut rng: rand::rngs::ThreadRng = thread_rng();
    mu[0] = rng.gen_range(0.0..2.0 * PI);

    let mut ci: f64 = 0.0;
    //オイラー法
    for k in 0..time.simulation_time - 1 {
        //濃度の更新
        c_vec.pop_front();
        c_vec.push_back(concentration(&constant, r[0][0], r[0][1]));

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

        //走化性能指数の計算
        ci += ((r[0][0] - constant.x_peak).powi(2) + (r[0][1] - constant.y_peak).powi(2)).sqrt();

        //更新
        for i in 0..8 {
            y[0][i] = y[1][i];
        }
        mu[0] = mu[1];
        for i in 0..2 {
            r[0][i] = r[1][i];
        }
    }
    //走化性能指数の計算
    ci += ((r[0][0] - constant.x_peak).powi(2) + (r[0][1] - constant.y_peak).powi(2)).sqrt();
    ci = 1.0
        - ci / (constant.x_peak.powi(2) + constant.y_peak.powi(2)).sqrt()
            / constant.simulation_time
            * constant.dt;
    if ci < 0.0 {
        ci = 0.0
    }
    ci
}

// chemotaxis_index関数は化学的走化性指数を計算します(線虫の行動に近くないと減点する)
pub fn chemotaxis_index_wave_check(gene: &Gene, setting: &Setting) -> f64 {
    //定数
    let constant: Const = setting.const_new();
    //遺伝子の受け渡し
    let weight: GeneConst = gene.scaling();
    //時間に関する定数をステップ数に変換
    let time: Time = time_new(&weight, &constant);

    //配列の宣言
    let mut y: [[f64; 8]; 2] = [[0.0; 8]; 2];
    let mut mu: [f64; 2] = [0.0; 2];
    let mut phi: f64;
    let mut r: [[f64; 2]; 2] = [[0.0; 2]; 2];
    let mut wave_check: f64 = 0.0;
    let mut wave_point: f64 = 0.0;

    //初期濃度の履歴生成
    let mut c_vec: VecDeque<f64> = VecDeque::new();
    for _ in 0..time.n_time + time.m_time {
        c_vec.push_back(concentration(&constant, 0.0, 0.0));
    }
    //運動ニューロンの初期活性を0～1の範囲でランダム化
    let _ = thread_rng().try_fill(&mut y[0][4..]);

    //ランダムな向きで配置
    let mut rng: rand::rngs::ThreadRng = thread_rng();
    mu[0] = rng.gen_range(0.0..2.0 * PI);

    let mut ci: f64 = 0.0;
    //オイラー法
    for k in 0..time.simulation_time - 1 {
        //濃度の更新
        c_vec.pop_front();
        c_vec.push_back(concentration(&constant, r[0][0], r[0][1]));

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
        phi = weight.w_nmj * (d - v);
        mu[1] = mu[0] + phi * constant.dt;

        //ピルエット
        if k % time.f_inv_time == time.f_inv_time - 1 {
            let mut rng: rand::rngs::ThreadRng = thread_rng();
            mu[1] = rng.gen_range(0.0..2.0 * PI);
        }

        //波打ちながら進んでいるかのチェック
        if k % time.periodic_time == time.periodic_time / 4 {
            wave_check = phi;
        } else if k % time.periodic_time == time.periodic_time * 3 / 4 && wave_check * phi > 0.0 {
            wave_point += 0.008;
        }

        //位置の更新
        r[1][0] = r[0][0] + constant.velocity * (mu[0]).cos() * constant.dt;
        r[1][1] = r[0][1] + constant.velocity * (mu[0]).sin() * constant.dt;

        //走化性能指数の計算
        ci += ((r[0][0] - constant.x_peak).powi(2) + (r[0][1] - constant.y_peak).powi(2)).sqrt();

        //更新
        for i in 0..8 {
            y[0][i] = y[1][i];
        }
        mu[0] = mu[1];
        for i in 0..2 {
            r[0][i] = r[1][i];
        }
    }
    //走化性能指数の計算
    ci += ((r[0][0] - constant.x_peak).powi(2) + (r[0][1] - constant.y_peak).powi(2)).sqrt();
    ci = 1.0
        - ci / (constant.x_peak.powi(2) + constant.y_peak.powi(2)).sqrt()
            / constant.simulation_time
            * constant.dt
        - wave_point;
    if ci < 0.0 {
        ci = 0.0
    }
    ci
}