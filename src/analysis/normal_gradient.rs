use crate::analysis::klinotaxis::klinotaxis_analysis;
use crate::analysis::setting::*;
use crate::genetic_algorithm::*;
use crate::simulation::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;

pub fn normal_gradient(
    r_vec: &[[f64; 2]],
    constant: &Const,
    time: &Time,
    periodic_number: &usize,
    delta: f64,
    mode: usize,
) -> Vec<f64> {
    // concentration関数の選択
    let concentration_fn: fn(&Const, f64, f64) -> f64 = match mode {
        0 => concentration,
        1 => gauss_concentration,
        2 => two_gauss_concentration,
        _ => panic!("Invalid mode: {}", mode),
    };

    let mut bearing_point: Vec<[[f64; 2]; 2]> = Vec::new();
    for i in 0..time.simulation_time - 2 * periodic_number * time.periodic_time {
        bearing_point.push([r_vec[i], r_vec[i + periodic_number * time.periodic_time]]);
    }
    let mut normal_gradient: Vec<f64> = Vec::new();

    // 法線ベクトル
    for bearing_point_item in &bearing_point {
        let bearing_vec: [f64; 2] = [
            -bearing_point_item[1][1] + bearing_point_item[0][1],
            bearing_point_item[1][0] - bearing_point_item[0][0],
        ];

        let magnitude_bearing: f64 = (bearing_vec[0].powi(2) + bearing_vec[1].powi(2)).sqrt();

        let normal_gradient_delta: [f64; 2] = [
            bearing_vec[0] / magnitude_bearing * delta,
            bearing_vec[1] / magnitude_bearing * delta,
        ];

        let x = bearing_point_item[0][0];
        let y = bearing_point_item[0][1];

        let normal = (concentration_fn(
            constant,
            x + normal_gradient_delta[0],
            y + normal_gradient_delta[1],
        ) - concentration_fn(constant, x, y))
            / delta;

        normal_gradient.push(normal);
    }

    normal_gradient
}

pub fn histgram_count_normal_gradient(
    nomal_gradient: &[f64],
    curving_rate: &[f64],
    bin_number: usize,
    concentration_gradient_max: f64,
) -> Vec<f64> {
    //ヒストグラムの種
    let mut curving_rate_mean: Vec<f64> = Vec::new();
    let step_size: f64 = 2.0 * concentration_gradient_max / bin_number as f64;

    for i in 0..bin_number {
        let mut mean: Vec<f64> = Vec::new();
        for (j, &nomal_gradient_value) in nomal_gradient.iter().enumerate() {
            if ((-concentration_gradient_max + (i as f64) * step_size) < nomal_gradient_value)
                && (nomal_gradient_value
                    < (-concentration_gradient_max + ((i + 1) as f64) * step_size))
            {
                mean.push(curving_rate[j]);
            }
        }
        // 平均値の計算
        let mean_value: f64 = if !mean.is_empty() {
            mean.iter().sum::<f64>() / mean.len() as f64
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };
        curving_rate_mean.push(mean_value);
    }

    curving_rate_mean
}

pub fn analysis_klinotaxis_nomal_gradient_errbar_std_max_min(
    result_ga: &[Ga],
    file_name: &Filename,
    liner_setting: &Setting,
    analysis_setting: &Analysissetting,
    analysis_use_function: &i32,
) {
    let hist_mean: Vec<Vec<f64>> = (0..analysis_setting.analysis_loop)
        .into_par_iter()
        .map(|_| {
            let result: Vec<Vec<f64>> = klinotaxis_analysis(
                &result_ga[analysis_setting.gene_number].gene,
                liner_setting,
                analysis_setting.mode,
                analysis_setting.periodic_number,
                analysis_setting.periodic_number_drain,
                analysis_setting.delta,
                analysis_use_function,
            );

            histgram_count_normal_gradient(
                &result[0],
                &result[1],
                analysis_setting.bin_number,
                analysis_setting.concentration_gradient_max,
            )
        })
        .collect::<Vec<Vec<f64>>>();

    // ヒストグラムの作成
    let mut normal_gradient_hist: Vec<f64> = Vec::new();
    let mut curving_rate_hist: Vec<f64> = Vec::new();
    let mut error_bar_std: Vec<f64> = Vec::new();
    let mut error_bar_max: Vec<f64> = Vec::new();
    let mut error_bar_min: Vec<f64> = Vec::new();

    let step_size: f64 =
        2.0 * analysis_setting.concentration_gradient_max / analysis_setting.bin_number as f64;

    for i in 0..analysis_setting.bin_number {
        let mut mean: Vec<f64> = Vec::new();
        for row in &hist_mean {
            if !row[i].is_nan() {
                mean.push(row[i]);
            }
        }

        // 平均値の計算
        let mean_value: f64 = if !mean.is_empty() {
            mean.iter().sum::<f64>() / mean.len() as f64
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        // エラーバーの計算
        let max: f64 = if !mean.is_empty() {
            mean.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        } else {
            f64::NAN
        };

        let min: f64 = if !mean.is_empty() {
            mean.iter().cloned().fold(f64::INFINITY, f64::min)
        } else {
            f64::NAN
        };

        // 標準偏差の計算
        let std: f64 = if !mean.is_empty() {
            let variance: f64 =
                mean.iter().map(|&x| (x - mean_value).powi(2)).sum::<f64>() / mean.len() as f64;
            variance.sqrt()
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        normal_gradient_hist
            .push(-analysis_setting.concentration_gradient_max + (i as f64) * step_size);
        curving_rate_hist.push(mean_value);
        error_bar_std.push(std);
        error_bar_max.push(max);
        error_bar_min.push(min);
    }

    // Open a file for writing
    let mut file: File = File::create(&file_name.nomal_gradient_vs_curving_rate_output).unwrap();

    // Iterate over the vectors and write each triplet of values to a line in the file
    for ((((nomal, curving_rate), std), max), min) in normal_gradient_hist
        .iter()
        .zip(curving_rate_hist.iter())
        .zip(error_bar_std.iter())
        .zip(error_bar_max.iter())
        .zip(error_bar_min.iter())
    {
        writeln!(
            file,
            "{}, {}, {}, {}, {}",
            nomal, curving_rate, std, max, min
        )
        .unwrap();
    }
}
