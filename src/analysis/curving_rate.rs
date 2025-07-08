use crate::simulation::*;

pub fn curving_rate_bear(r_vec: &[[f64; 2]], time: &Time, periodic_number: &usize) -> Vec<f64> {
    //  カービングレート計算用の座標リストを初期化
    let mut curving_point: Vec<[[f64; 2]; 2]> = Vec::new();
    for i in 0..time.simulation_time - periodic_number * time.periodic_time {
        curving_point.push([r_vec[i], r_vec[i + periodic_number * time.periodic_time]]);
    }

    //  ベクトル・距離リストを初期化
    let mut curving_vec: Vec<[[f64; 2]; 2]> = Vec::new();
    let mut curving_distance: Vec<f64> = Vec::new();

    for i in 0..time.simulation_time - 2 * periodic_number * time.periodic_time {
        let vec1 = [
            curving_point[i][1][0] - curving_point[i][0][0],
            curving_point[i][1][1] - curving_point[i][0][1],
        ];
        let vec2 = [
            curving_point[i + periodic_number * time.periodic_time][1][0]
                - curving_point[i + periodic_number * time.periodic_time][0][0],
            curving_point[i + periodic_number * time.periodic_time][1][1]
                - curving_point[i + periodic_number * time.periodic_time][0][1],
        ];
        curving_vec.push([vec1, vec2]);

        let distance = vec1.iter().map(|x| x.powi(2)).sum::<f64>().sqrt()
            + vec2.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();
        curving_distance.push(distance);
    }

    //  カービング角度（角度差）リストを初期化
    let mut curving_angle: Vec<f64> = Vec::new();

    for vec_pair in &curving_vec {
        let dot = vec_pair[0][0] * vec_pair[1][0] + vec_pair[0][1] * vec_pair[1][1];
        let mag1 = (vec_pair[0][0].powi(2) + vec_pair[0][1].powi(2)).sqrt();
        let mag2 = (vec_pair[1][0].powi(2) + vec_pair[1][1].powi(2)).sqrt();
        let angle_rad = (dot / (mag1 * mag2)).acos();
        let mut angle_deg = angle_rad.to_degrees();

        let cross = vec_pair[0][0] * vec_pair[1][1] - vec_pair[0][1] * vec_pair[1][0];
        if cross < 0.0 {
            angle_deg *= -1.0;
        }
        curving_angle.push(angle_deg);
    }

    //  カービングレートの計算
    curving_angle
        .iter()
        .zip(&curving_distance)
        .map(|(angle, dist)| angle / dist)
        .collect()
}
