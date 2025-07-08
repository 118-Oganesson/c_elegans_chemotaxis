// 遺伝子のデータを表現するための構造体Gene
#[derive(Clone, Debug)]
pub struct Gene {
    pub gene: Vec<f64>, // 遺伝子データを格納するベクタ
}

// 遺伝子の定数を格納するための構造体GeneConst
pub struct GeneConst {
    pub n: f64,           // 感覚ニューロンの時間スケール
    pub m: f64,           // 運動ニューロンの時間スケール
    pub theta: [f64; 8],  // ニューロンの閾値
    pub w_on: [f64; 8],   // 感覚ニューロンONのシナプス結合重み
    pub w_off: [f64; 8],  // 感覚ニューロンOFFのシナプス結合重み
    pub w: [[f64; 8]; 8], // 介在ニューロンと運動ニューロンのシナプス結合重み
    pub g: [[f64; 8]; 8], // 介在ニューロンと運動ニューロンのギャップ結合重み
    pub w_osc: [f64; 8],  // 運動ニューロンに入る振動成分の重み
    pub w_nmj: f64,       // 回転角度の重み
}

impl Gene {
    // Gene構造体のデータを基にGeneConst構造体を生成するメソッド
    pub fn scaling(&self) -> GeneConst {
        // 内部関数range: 遺伝子データを指定範囲にスケーリングする
        fn range(gene: &f64, min: f64, max: f64) -> f64 {
            (gene + 1.0) / 2.0 * (max - min) + min
        }

        // ニューロンの閾値 [-15, 15]の範囲でスケーリング
        let mut theta: [f64; 8] = [0.0; 8];
        theta[0] = range(&self.gene[2], -15.0, 15.0);
        theta[1] = range(&self.gene[3], -15.0, 15.0);
        theta[2] = range(&self.gene[4], -15.0, 15.0);
        theta[3] = range(&self.gene[5], -15.0, 15.0);
        theta[4] = range(&self.gene[6], -15.0, 15.0);
        theta[5] = theta[4];
        theta[6] = range(&self.gene[7], -15.0, 15.0);
        theta[7] = theta[6];

        // 感覚ニューロンONのシナプス結合重み [-15, 15]の範囲でスケーリング
        let mut w_on: [f64; 8] = [0.0; 8];
        w_on[0] = range(&self.gene[8], -15.0, 15.0);
        w_on[1] = range(&self.gene[9], -15.0, 15.0);

        // 感覚ニューロンOFFのシナプス結合重み [-15, 15]の範囲でスケーリング
        let mut w_off: [f64; 8] = [0.0; 8];
        w_off[0] = range(&self.gene[10], -15.0, 15.0);
        w_off[1] = range(&self.gene[11], -15.0, 15.0);

        // 介在ニューロンと運動ニューロンのシナプス結合重み [-15, 15]の範囲でスケーリング
        let mut w: [[f64; 8]; 8] = [[0.0; 8]; 8];
        w[0][2] = range(&self.gene[12], -15.0, 15.0);
        w[1][3] = range(&self.gene[13], -15.0, 15.0);
        w[2][4] = range(&self.gene[14], -15.0, 15.0);
        w[2][5] = w[2][4];
        w[3][6] = range(&self.gene[15], -15.0, 15.0);
        w[3][7] = w[3][6];
        w[4][4] = range(&self.gene[16], -15.0, 15.0);
        w[5][5] = w[4][4];
        w[6][6] = range(&self.gene[17], -15.0, 15.0);
        w[7][7] = w[6][6];

        // 介在ニューロンと運動ニューロンのギャップ結合重み [0, 2.5]の範囲でスケーリング
        let mut g: [[f64; 8]; 8] = [[0.0; 8]; 8];
        g[0][1] = range(&self.gene[18], 0.0, 2.5);
        g[1][0] = g[0][1];
        g[2][3] = range(&self.gene[19], 0.0, 2.5);
        g[3][2] = g[2][3];

        // 運動ニューロンに入る振動成分の重み [0, 15]の範囲でスケーリング
        let mut w_osc: [f64; 8] = [0.0; 8];
        w_osc[4] = range(&self.gene[20], 0.0, 15.0);
        w_osc[7] = w_osc[4];
        w_osc[5] = -w_osc[4];
        w_osc[6] = -w_osc[4];

        // GeneConst構造体を生成して返す
        GeneConst {
            n: range(&self.gene[0], 0.1, 4.2), // 感覚ニューロン時間スケールのスケーリング [0.1, 4.2]
            m: range(&self.gene[1], 0.1, 4.2), // 運動ニューロン時間スケーリング [0.1, 4.2]
            theta,
            w_on,
            w_off,
            w,
            g,
            w_osc,
            w_nmj: range(&self.gene[21], 1.0, 3.0), // 回転角度の重み [1, 3]
        }
    }
}

// Gene構造体をf64のイテレータから生成するためのFromIteratorトレイトの実装
impl FromIterator<f64> for Gene {
    fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
        let gene: Vec<f64> = iter.into_iter().collect(); // イテレータからベクタを生成
        Gene { gene } // Gene構造体を生成して返す
    }
}
