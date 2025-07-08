use crate::analysis::setting::*;
use crate::genetic_algorithm::*;
use crate::simulation::*;
use toml::Value;

pub fn read_result() -> Vec<Ga> {
    // klinotaxis_analysis.toml ファイルを読み込む
    let toml_str: String =
        std::fs::read_to_string("klinotaxis_analysis_setting.toml").expect("Failed to read file");
    let value: Value = toml::from_str(&toml_str).expect("Failed to parse TOML");
    let file_name: Filename = value["file_name"]
        .clone()
        .try_into()
        .expect("Failed to parse Setting");

    let result_json: String = std::fs::read_to_string(file_name.read_result_json)
        .expect("ファイルの読み込みに失敗しました");

    // JSON 文字列を Vec<Gajson> にデシリアライズ
    let result_gajson: Vec<Gajson> =
        serde_json::from_str(&result_json).expect("JSONのデシリアライズに失敗しました");

    // Gajson を Ga に変換
    let result: Vec<Ga> = result_gajson
        .into_iter()
        .map(|gajson: Gajson| Ga {
            value: gajson.value,
            gene: Gene { gene: gajson.gene },
        })
        .collect();
    result
}
