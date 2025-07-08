use serde::{Deserialize, Serialize};

// ファイル名関連の情報を格納する構造体Filename
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Filename {
    pub read_result_json: String, // 解析結果の JSON ファイルのパス
    pub bearing_vs_curving_rate_output: String, // Bearing vs Turning Bias の出力ファイルのパス
    pub nomal_gradient_vs_curving_rate_output: String, // Nomal Gradient vs Turning Bias の出力ファイルのパス
    pub translational_gradient_vs_curving_rate_output: String, // Translational Gradient vs Turning Bias の出力ファイルのパス
}

// 解析に使用される関数の情報を格納する構造体Analysisusefunction
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Analysisusefunction {
    pub mode: Vec<i32>, // 解析関数のモードを示す番号のリスト
}

// 解析に使用される遺伝子の情報を格納する構造体Analysisusegene
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Analysisusegene {
    pub mode: usize,                   // 解析に使用する遺伝子の設定
    pub gene_numbers: Vec<usize>,      // 解析する個別の遺伝子の番号のリスト
    pub gene_number_range: Vec<usize>, // 解析する遺伝子の範囲を示すリスト（開始番号と終了番号を含む）
}

// 解析設定の情報を格納する構造体Analysissetting
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Analysissetting {
    pub gene_number: usize,              // 解析対象の遺伝子の番号
    pub mode: usize, // 濃度関数の選択（mode=0:liner, mode=1:gauss, mode=2:two_gauss）
    pub analysis_loop: usize, // 統計を取る回数
    pub periodic_number: usize, // 解析に使用するデータの周期数
    pub periodic_number_drain: usize, // 解析において無視する初期の周期数
    pub bin_range: usize, // ヒストグラムのビンの幅（bearing）
    pub delta: f64,  // 濃度勾配の2点間距離
    pub bin_number: usize, // ヒストグラムのビンの数（濃度勾配）
    pub concentration_gradient_max: f64, // 濃度勾配の最大値
}
