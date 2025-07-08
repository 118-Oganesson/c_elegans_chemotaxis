use crate::simulation::gene::*;
use crate::simulation::setting::*;
use crate::simulation::time::*;
use std::collections::VecDeque;
use std::f64::consts::PI;

// concentration関数は与えられた位置(x, y)における濃度を計算します
#[inline]
pub fn concentration(constant: &Const, x: f64, y: f64) -> f64 {
    constant.alpha * ((x - constant.x_peak).powi(2) + (y - constant.y_peak).powi(2)).sqrt()
}

// gauss_concentration関数はガウス関数に基づいて与えられた位置(x, y)における濃度を計算します
#[inline]
pub fn gauss_concentration(constant: &Const, x: f64, y: f64) -> f64 {
    constant.c_0
        * (-((x - constant.x_peak).powi(2) + (y - constant.y_peak).powi(2))
            / (2.0 * constant.lambda.powi(2)))
        .exp()
}

// two_gauss_concentration関数は2つのガウス関数に基づいて与えられた位置(x, y)における濃度を計算します
#[inline]
pub fn two_gauss_concentration(constant: &Const, x: f64, y: f64) -> f64 {
    constant.c_0
        * ((-((x - constant.x_peak).powi(2) + (y - constant.y_peak).powi(2))
            / (2.0 * constant.lambda.powi(2)))
        .exp()
            - (-((x + constant.x_peak).powi(2) + (y + constant.y_peak).powi(2))
                / (2.0 * constant.lambda.powi(2)))
            .exp())
}

// sigmoid関数はシグモイド関数を計算します
#[inline]
pub fn sigmoid(x: f64) -> f64 {
    if x > 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

// y_osc関数は振動成分y_oscを計算します
#[inline]
pub fn y_osc(constant: &Const, time: f64) -> f64 {
    (2.0 * PI * time / constant.periodic_time).sin()
}

// y_on_off関数は濃度履歴からON細胞とOFF細胞の出力y_on、y_offを計算します
#[inline]
pub fn y_on_off(
    constant: &Const,
    weight: &GeneConst,
    time: &Time,
    c_vec: &VecDeque<f64>,
) -> [f64; 2] {
    let y_on: f64 = c_vec
        .iter()
        .enumerate()
        .map(|(i, &value)| {
            if i < time.m_time {
                -value / weight.m
            } else {
                value / weight.n
            }
        })
        .sum();

    if y_on < 0.0 {
        [0.0, -y_on * 100.0 * constant.dt]
    } else {
        [y_on * 100.0 * constant.dt, 0.0]
    }
}
