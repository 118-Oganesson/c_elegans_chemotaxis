use crate::analysis::klinotaxis::klinotaxis_analysis;
use crate::analysis::setting::*;
use crate::genetic_algorithm::*;
use crate::simulation::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;

pub fn translational_gradient(
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
    let mut translational_gradient: Vec<f64> = Vec::new();

    // 並進方向のベクトル
    for bearing_point_item in &bearing_point {
        let bearing_vec: [f64; 2] = [
            bearing_point_item[1][0] - bearing_point_item[0][0],
            bearing_point_item[1][1] - bearing_point_item[0][1],
        ];

        let magnitude_bearing: f64 = (bearing_vec[0].powi(2) + bearing_vec[1].powi(2)).sqrt();

        let translational_gradient_delta: [f64; 2] = [
            bearing_vec[0] / magnitude_bearing * delta,
            bearing_vec[1] / magnitude_bearing * delta,
        ];

        let x = bearing_point_item[0][0];
        let y = bearing_point_item[0][1];

        let translational = (concentration_fn(
            constant,
            x + translational_gradient_delta[0],
            y + translational_gradient_delta[1],
        ) - concentration_fn(constant, x, y))
            / delta;

        translational_gradient.push(translational);
    }

    translational_gradient
}

pub fn histgram_count_translational_gradient(
    translational_gradient: &[f64],
    curving_rate: &[f64],
    bin_number: usize,
    concentration_gradient_max: f64,
) -> Vec<Vec<f64>> {
    // ヒストグラムの種
    let mut curving_rate_mean: Vec<f64> = Vec::new();
    let mut positive_mean: Vec<f64> = Vec::new();
    let mut negative_mean: Vec<f64> = Vec::new();
    let step_size: f64 = 2.0 * concentration_gradient_max / bin_number as f64;

    for i in 0..bin_number {
        let mut mean: Vec<f64> = Vec::new();

        for (j, &translational_gradient_value) in translational_gradient.iter().enumerate() {
            if ((-concentration_gradient_max + (i as f64) * step_size)
                < translational_gradient_value)
                && (translational_gradient_value
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

        // 正の値のみの平均
        let positive_values: Vec<f64> = mean.iter().cloned().filter(|&x| x > 0.0).collect();
        let positive_mean_value: f64 = if !positive_values.is_empty() {
            positive_values.iter().sum::<f64>() / positive_values.len() as f64
        } else {
            f64::NAN
        };
        positive_mean.push(positive_mean_value);

        // 負の値のみの平均
        let negative_values: Vec<f64> = mean.iter().cloned().filter(|&x| x < 0.0).collect();
        let negative_mean_value: f64 = if !negative_values.is_empty() {
            negative_values.iter().sum::<f64>() / negative_values.len() as f64
        } else {
            f64::NAN
        };
        negative_mean.push(negative_mean_value);
    }

    // 結果を返す
    let result: Vec<Vec<f64>> = vec![curving_rate_mean, positive_mean, negative_mean];

    result
}

pub fn analysis_klinotaxis_translational_gradient_errbar_std_positive_negative(
    result_ga: &[Ga],
    file_name: &Filename,
    liner_setting: &Setting,
    analysis_setting: &Analysissetting,
    analysis_use_function: &i32,
) {
    let hist_mean: Vec<Vec<Vec<f64>>> = (0..analysis_setting.analysis_loop)
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

            histgram_count_translational_gradient(
                &result[0],
                &result[1],
                analysis_setting.bin_number,
                analysis_setting.concentration_gradient_max,
            )
        })
        .collect::<Vec<Vec<Vec<f64>>>>();

    // ヒストグラムの作成
    let mut translational_gradient_hist: Vec<f64> = Vec::new();
    let mut curving_rate_hist: Vec<f64> = Vec::new();
    let mut curving_rate_positive_hist: Vec<f64> = Vec::new();
    let mut curving_rate_negative_hist: Vec<f64> = Vec::new();
    let mut error_bar_std: Vec<f64> = Vec::new();
    let mut error_bar_positive_std: Vec<f64> = Vec::new();
    let mut error_bar_negative_std: Vec<f64> = Vec::new();

    let step_size: f64 =
        2.0 * analysis_setting.concentration_gradient_max / analysis_setting.bin_number as f64;

    for i in 0..analysis_setting.bin_number {
        let mut mean: Vec<f64> = Vec::new();
        let mut positive_mean: Vec<f64> = Vec::new();
        let mut negative_mean: Vec<f64> = Vec::new();

        for row in &hist_mean {
            if !row[0][i].is_nan() {
                mean.push(row[0][i]);
            }
        }

        for row in &hist_mean {
            if !row[1][i].is_nan() {
                positive_mean.push(row[1][i]);
            }
        }

        for row in &hist_mean {
            if !row[2][i].is_nan() {
                negative_mean.push(row[2][i]);
            }
        }

        // 平均値の計算
        let mean_value: f64 = if !mean.is_empty() {
            mean.iter().sum::<f64>() / mean.len() as f64
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        let positive_mean_value: f64 = if !positive_mean.is_empty() {
            positive_mean.iter().sum::<f64>() / positive_mean.len() as f64
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        let negative_mean_value: f64 = if !negative_mean.is_empty() {
            negative_mean.iter().sum::<f64>() / negative_mean.len() as f64
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        // 標準偏差の計算
        let std: f64 = if !mean.is_empty() {
            let variance: f64 =
                mean.iter().map(|&x| (x - mean_value).powi(2)).sum::<f64>() / mean.len() as f64;
            variance.sqrt()
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        let positive_std: f64 = if !positive_mean.is_empty() {
            let variance: f64 = positive_mean
                .iter()
                .map(|&x| (x - positive_mean_value).powi(2))
                .sum::<f64>()
                / positive_mean.len() as f64;
            variance.sqrt()
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        let negative_std: f64 = if !negative_mean.is_empty() {
            let variance: f64 = negative_mean
                .iter()
                .map(|&x| (x - negative_mean_value).powi(2))
                .sum::<f64>()
                / negative_mean.len() as f64;
            variance.sqrt()
        } else {
            f64::NAN // もし mean が空なら NaN をセット
        };

        translational_gradient_hist
            .push(-analysis_setting.concentration_gradient_max + (i as f64) * step_size);
        curving_rate_hist.push(mean_value);
        curving_rate_positive_hist.push(positive_mean_value);
        curving_rate_negative_hist.push(negative_mean_value);
        error_bar_std.push(std);
        error_bar_positive_std.push(positive_std);
        error_bar_negative_std.push(negative_std);
    }

    // Open a file for writing
    let mut file: File =
        File::create(&file_name.translational_gradient_vs_curving_rate_output).unwrap();

    // Iterate over the vectors and write each triplet of values to a line in the file
    for (
        (
            ((((translational, curving_rate), std), positive_curving_rate), positive_std),
            negative_curving_rate,
        ),
        negative_std,
    ) in translational_gradient_hist
        .iter()
        .zip(curving_rate_hist.iter())
        .zip(error_bar_std.iter())
        .zip(curving_rate_positive_hist.iter())
        .zip(error_bar_positive_std.iter())
        .zip(curving_rate_negative_hist.iter())
        .zip(error_bar_negative_std.iter())
    {
        writeln!(
            file,
            "{}, {}, {}, {}, {}, {}, {}",
            translational,
            curving_rate,
            std,
            positive_curving_rate,
            positive_std,
            negative_curving_rate,
            negative_std,
        )
        .unwrap();
    }
}
