use crate::analysis::bearing::*;
use crate::analysis::normal_gradient::*;
use crate::analysis::read_result::*;
use crate::analysis::setting::*;
use crate::analysis::translational_gradient::*;
use crate::genetic_algorithm::*;
use crate::simulation::*;
use std::fs;
use std::path::Path;
use toml::Value;

pub fn analytical_control() {
    // 結果読み込み
    let result_ga: Vec<Ga> = read_result();

    // TOMLファイル読み込みと解析
    let toml_str = match fs::read_to_string("klinotaxis_analysis_setting.toml") {
        Ok(content) => content,
        Err(e) => {
            eprintln!("Error reading klinotaxis_analysis_setting.toml: {}", e);
            return;
        }
    };

    let value: Value = match toml::from_str(&toml_str) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("Error parsing TOML: {}", e);
            return;
        }
    };

    // 設定の読み出し（各フィールドで丁寧にエラー処理）
    let file_name_base: Filename = value
        .get("file_name")
        .and_then(|v| v.clone().try_into().ok())
        .unwrap_or_else(|| {
            panic!("Missing or invalid [file_name] section in TOML");
        });

    let liner_setting: Setting = value
        .get("liner_setting")
        .and_then(|v| v.clone().try_into().ok())
        .unwrap_or_else(|| {
            panic!("Missing or invalid [liner_setting] section in TOML");
        });

    let mut analysis_setting: Analysissetting = value
        .get("analysis_setting")
        .and_then(|v| v.clone().try_into().ok())
        .unwrap_or_else(|| {
            panic!("Missing or invalid [analysis_setting] section in TOML");
        });

    let analysis_use_gene: Analysisusegene = value
        .get("analysis_use_gene")
        .and_then(|v| v.clone().try_into().ok())
        .unwrap_or_else(|| {
            panic!("Missing or invalid [analysis_use_gene] section in TOML");
        });

    let analysis_use_function: Analysisusefunction = value
        .get("analysis_use_function")
        .and_then(|v| v.clone().try_into().ok())
        .unwrap_or_else(|| {
            panic!("Missing or invalid [analysis_use_function] section in TOML");
        });

    // 対象遺伝子リスト作成
    let gene_range: Vec<usize> = match analysis_use_gene.mode {
        0 => analysis_use_gene.gene_numbers.clone(),
        1 => (analysis_use_gene.gene_number_range[0]..=analysis_use_gene.gene_number_range[1])
            .collect(),
        _ => {
            eprintln!("Invalid gene selection mode: {}", analysis_use_gene.mode);
            return;
        }
    };

    // 各遺伝子ごとの分析処理
    for gene_index in gene_range {
        analysis_setting.gene_number = gene_index;
        let file_name = generate_file_names(&file_name_base, gene_index);

        // 出力先ディレクトリの作成（必要なら）
        ensure_parent_dirs(&file_name);

        for mode in &analysis_use_function.mode {
            println!("Analyzing gene {} with mode {}", gene_index, mode);
            match mode {
                0 => analysis_klinotaxis_bearing_errbar_std_max_min(
                    &result_ga,
                    &file_name,
                    &liner_setting,
                    &analysis_setting,
                    mode,
                ),
                1 => analysis_klinotaxis_nomal_gradient_errbar_std_max_min(
                    &result_ga,
                    &file_name,
                    &liner_setting,
                    &analysis_setting,
                    mode,
                ),
                2 => analysis_klinotaxis_translational_gradient_errbar_std_positive_negative(
                    &result_ga,
                    &file_name,
                    &liner_setting,
                    &analysis_setting,
                    mode,
                ),
                _ => {
                    eprintln!("Unknown analysis mode: {}", mode);
                }
            }
        }
    }
}

/// ファイル名構造体を遺伝子番号付きで生成
fn generate_file_names(base: &Filename, index: usize) -> Filename {
    Filename {
        bearing_vs_curving_rate_output: format!(
            "{}_{}.txt",
            base.bearing_vs_curving_rate_output, index
        ),
        nomal_gradient_vs_curving_rate_output: format!(
            "{}_{}.txt",
            base.nomal_gradient_vs_curving_rate_output, index
        ),
        translational_gradient_vs_curving_rate_output: format!(
            "{}_{}.txt",
            base.translational_gradient_vs_curving_rate_output, index
        ),
        ..base.clone()
    }
}

/// ファイルパスの親ディレクトリが存在しなければ作成
fn ensure_parent_dirs(file_name: &Filename) {
    let paths = vec![
        &file_name.bearing_vs_curving_rate_output,
        &file_name.nomal_gradient_vs_curving_rate_output,
        &file_name.translational_gradient_vs_curving_rate_output,
    ];

    for path_str in paths {
        if let Some(parent) = Path::new(path_str).parent() {
            if !parent.exists() {
                if let Err(e) = fs::create_dir_all(parent) {
                    eprintln!("Failed to create directory {:?}: {}", parent, e);
                }
            }
        }
    }
}
