use crate::genetic_algorithm::crossover::*;
use crate::genetic_algorithm::evaluate::*;
use crate::genetic_algorithm::ga::*;
use crate::genetic_algorithm::mutation::*;
use crate::genetic_algorithm::population::*;
use crate::genetic_algorithm::setting::*;
use crate::simulation::*;
use rand::{rng, Rng};
use std::fs::File;
use std::io::prelude::*;
use toml::Value;

// genetic_algorithm_constrain_aiy_aiz関数は生物学的な制約を加えた遺伝的アルゴリズムを用いて線虫の生理学的なパラメータを最適化します
#[allow(dead_code)]
pub fn genetic_algorithm_constrain_aiy_aiz() {
    // tomlファイルを読み込む
    let toml_str: String =
        std::fs::read_to_string("genetic_algorithm_setting.toml").expect("Failed to read file");
    let value: Value = toml::from_str(&toml_str).expect("Failed to parse TOML");

    // 各セクションを取り出す
    let ga_setting: Gasetting = value["Ga_setting"]
        .clone()
        .try_into()
        .expect("Failed to parse Gasetting");
    let setting: Setting = value["setting"]
        .clone()
        .try_into()
        .expect("Failed to parse Setting");
    let testing: Setting = value["testing"]
        .clone()
        .try_into()
        .expect("Failed to parse Testing");

    // 途中結果の書き出し用のtxtファイルの作成
    let mut file = File::create("Result.txt").unwrap();

    //遺伝的アルゴリズムの結果
    let mut result: Vec<Ga> = Vec::new();

    for count in 0..ga_setting.ga_count {
        //初期集団を生成
        let mut population: Vec<Ga> =
            population_new_aiy_aiz_negative(ga_setting.gen_size, ga_setting.pop_size);

        //個体の評価(version:0は通常、version:1は波打つかチェックしている)
        let mut evaluate: Vec<Ga> =
            evaluate_fitness(&mut population, &setting, ga_setting.average, 1);

        //個体をvalueの値によって降順でsort
        evaluate.sort_by(|a: &Ga, b: &Ga| b.value.partial_cmp(&a.value).unwrap());

        println!(
            "{:3}_Gen: {:03}, Fitness_1: {:.5}",
            count + 1,
            0,
            evaluate[0].value
        );
        println!("              Fitness_2: {:.5}", evaluate[1].value);
        println!("              Fitness_3: {:.5}", evaluate[2].value);
        println!();

        for i in 0..ga_setting.n_gen {
            //選択
            let select: Vec<Ga> = evaluate.iter().take(ga_setting.sel_top).cloned().collect();

            //交叉
            let mut mate: Vec<Ga> = Vec::new();
            let clone: Vec<Ga> = select.clone();
            let mut rng = rng();
            for i in (0..clone.len()).step_by(2) {
                if i + 1 < clone.len() && rng.random::<f64>() < ga_setting.mat_pb {
                    let parent1: &Ga = &clone[i];
                    let parent2: &Ga = &clone[i + 1];
                    let child: Vec<Ga> = two_point_crossover(parent1, parent2);
                    mate.extend(child);
                }
            }

            //変異
            let mut mutant: Vec<Ga> = Vec::new();
            let clone: Vec<Ga> = select.clone();
            for ind in clone.iter() {
                if rng.random::<f64>() < ga_setting.mut_pb {
                    mutant.push(mutation_aiy_aiz_negative(ind, 0.4, (0.0, 0.05)));
                }
            }

            if i % ga_setting.re_val == 0 {
                //再評価
                let mut offspring: Vec<Ga> = Vec::new();
                offspring.extend(select);
                offspring.extend(mate);
                offspring.extend(mutant);
                let population: Vec<Ga> = population_new_aiy_aiz_negative(
                    ga_setting.gen_size,
                    ga_setting.pop_size - offspring.len(),
                );
                offspring.extend(population);
                evaluate.clear();
                let offspring_evaluate: Vec<Ga> =
                    evaluate_fitness(&mut offspring, &setting, ga_setting.average, 1);
                evaluate.extend(offspring_evaluate);
                evaluate.sort_by(|a, b| b.value.partial_cmp(&a.value).unwrap());
            } else {
                //通常評価
                let mut offspring: Vec<Ga> = Vec::new();
                offspring.extend(mate);
                offspring.extend(mutant);
                let population: Vec<Ga> = population_new_aiy_aiz_negative(
                    ga_setting.gen_size,
                    ga_setting.pop_size - ga_setting.sel_top - offspring.len(),
                );
                offspring.extend(population);
                evaluate.clear();
                let offspring_evaluate: Vec<Ga> =
                    evaluate_fitness(&mut offspring, &setting, ga_setting.average, 1);
                evaluate.extend(offspring_evaluate);
                evaluate.extend(select);
                evaluate.sort_by(|a, b| b.value.partial_cmp(&a.value).unwrap());
            }

            println!(
                "{:3}_Gen: {:03}, Fitness_1: {:.5}",
                count + 1,
                i + 1,
                evaluate[0].value
            );
            println!("              Fitness_2: {:.5}", evaluate[1].value);
            println!("              Fitness_3: {:.5}", evaluate[2].value);
            println!();
        }

        //最も優秀な個体を結果に格納
        result.push(evaluate[0].clone());

        // 途中結果の書き出し
        let result_gajson: Gajson = Gajson {
            value: evaluate[0].value,
            gene: evaluate[0].gene.gene.clone(),
        };
        let json_string = serde_json::to_string_pretty(&result_gajson).unwrap();
        file.write_all(json_string.as_bytes()).unwrap();
        writeln!(file, ",").unwrap();
    }

    //正しいCIを用いて結果を評価する
    let mut result_evaluate: Vec<Ga> =
        evaluate_fitness(&mut result, &testing, ga_setting.average, 0);
    result_evaluate.sort_by(|a, b| b.value.partial_cmp(&a.value).unwrap());

    //Gajsonに変換
    let result_evaluate_gajson: Vec<Gajson> = result_evaluate
        .into_iter()
        .map(|ga: Ga| Gajson {
            value: ga.value,
            gene: ga.gene.gene,
        })
        .collect();

    //JSON文字列にシリアライズ
    let result_json: String = serde_json::to_string_pretty(&result_evaluate_gajson).unwrap();

    //JSON文字列をファイルに書き込む
    let mut file: File =
        File::create("Result_aiy_aiz_negative.json").expect("ファイルの作成に失敗しました");
    file.write_all(result_json.as_bytes())
        .expect("ファイルへの書き込みに失敗しました");
}
