use crate::analysis::klinotaxis::klinotaxis_analysis;
use crate::analysis::setting::*;
use crate::genetic_algorithm::*;
use crate::simulation::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;

pub fn bearing(
    r_vec: &[[f64; 2]],
    constant: &Const,
    time: &Time,
    periodic_number: &usize,
) -> Vec<f64> {
    let mut bearing_point: Vec<[[f64; 2]; 2]> = Vec::new();
    for i in 0..time.simulation_time - 2 * periodic_number * time.periodic_time {
        bearing_point.push([r_vec[i], r_vec[i + periodic_number * time.periodic_time]]);
    }

    let peak_vec: [f64; 2] = [constant.x_peak, constant.y_peak];
    let mut bearing: Vec<f64> = Vec::new();

    for bearing_point_item in &bearing_point {
        let bearing_vec: [f64; 2] = [
            bearing_point_item[1][0] - bearing_point_item[0][0],
            bearing_point_item[1][1] - bearing_point_item[0][1],
        ];
        let dot_product: f64 = peak_vec[0] * bearing_vec[0] + peak_vec[1] * bearing_vec[1];
        let magnitude_peak: f64 = (peak_vec[0].powi(2) + peak_vec[1].powi(2)).sqrt();
        let magnitude_bearing: f64 = (bearing_vec[0].powi(2) + bearing_vec[1].powi(2)).sqrt();
        let angle_radian: f64 = (dot_product / (magnitude_peak * magnitude_bearing)).acos();
        let mut angle_degrees: f64 = angle_radian.to_degrees();
        let cross_product: f64 = peak_vec[0] * bearing_vec[1] - peak_vec[1] * bearing_vec[0];
        if cross_product > 0.0 {
            angle_degrees *= -1.0;
        }
        bearing.push(angle_degrees);
    }

    bearing
}

pub fn histgram_count_bearing(bearing: &[f64], curving_rate: &[f64], bin_range: usize) -> Vec<f64> {
    //ヒストグラムの種
    let mut curving_rate_mean: Vec<f64> = Vec::new();

    for i in (-180..180).step_by(bin_range) {
        let mut mean: Vec<f64> = Vec::new();
        for (j, &bearing_value) in bearing.iter().enumerate() {
            if ((i as f64) < bearing_value) && (bearing_value < (i as f64 + bin_range as f64)) {
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

pub fn analysis_klinotaxis_bearing_errbar_std_max_min(
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

            histgram_count_bearing(&result[0], &result[1], analysis_setting.bin_range)
        })
        .collect::<Vec<Vec<f64>>>();

    // ヒストグラムの作成
    let mut bearing_hist: Vec<f64> = Vec::new();
    let mut curving_rate_hist: Vec<f64> = Vec::new();
    let mut error_bar_std: Vec<f64> = Vec::new();
    let mut error_bar_max: Vec<f64> = Vec::new();
    let mut error_bar_min: Vec<f64> = Vec::new();

    for (count, i) in (-180..180).step_by(analysis_setting.bin_range).enumerate() {
        let mut mean: Vec<f64> = Vec::new();
        for row in &hist_mean {
            if !row[count].is_nan() {
                mean.push(row[count]);
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

        bearing_hist.push(i as f64);
        curving_rate_hist.push(mean_value);
        error_bar_std.push(std);
        error_bar_max.push(max);
        error_bar_min.push(min);
    }

    // Open a file for writing
    let mut file: File = File::create(&file_name.bearing_vs_curving_rate_output).unwrap();

    // Iterate over the vectors and write each triplet of values to a line in the file
    for ((((bearing, curving_rate), std), max), min) in bearing_hist
        .iter()
        .zip(curving_rate_hist.iter())
        .zip(error_bar_std.iter())
        .zip(error_bar_max.iter())
        .zip(error_bar_min.iter())
    {
        writeln!(
            file,
            "{}, {}, {}, {}, {}",
            bearing, curving_rate, std, max, min
        )
        .unwrap();
    }
}
