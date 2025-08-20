# `simulation` ディレクトリの説明

このディレクトリには、遺伝子パラメータセットに基づいて行った *C. elegans* の行動シミュレーション結果が複数のテキストファイルとして格納されています。各ファイルには、さまざまな条件下で計算された行動指標（例：屈曲率）の統計解析が記録されています。ファイルは名前とサブディレクトリごとに整理されています。

---

## 📁 ファイル & ディレクトリ概要

### 1. `with_constraining_AIY-AIZ/`

* **内容**: AIY–AIZ シナプスに抑制制約を課した条件で遺伝子が最適化された、最良の個体に関するシミュレーション結果。
* **ファイル**:

  * `bearing_vs_curving_rate_0.txt`
  * `normal_gradient_vs_curving_rate_0.txt`
  * `translational_gradient_vs_curving_rate_0.txt`
* **カラム**:

  * **bearing\_vs\_curving\_rate**: `[方位, 屈曲率, 標準偏差, 最大値, 最小値]`
  * **normal\_gradient\_vs\_curving\_rate**: `[法線勾配, 屈曲率, 標準偏差, 最大値, 最小値]`
  * **translational\_gradient\_vs\_curving\_rate**: `[移動勾配, 屈曲率, 標準偏差, 最大値, 最小値, 旋回バイアス（+）, 標準偏差（+）, 旋回バイアス（–）, 標準偏差（–）]`

---

## 📁 `concentration_memory/`

塩濃度記憶を組み込んだ走化性シミュレーションの解析結果を格納。

### 2. `Result_aiy_aiz_negative_0/`

* **内容**: `Result_aiy_aiz_negative.json` の遺伝子（インデックス 0）の解析データ。ASER–AIY シナプスを段階的に変化させたもの。
* **ファイル**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` ～ `n_vs_c_15.txt`

---

## 📁 `starvation/`

飢餓条件下でのシミュレーション解析を格納。

### 🔸 `synapse/`

#### 3. `Result_aiz_smb_0/`

* **内容**: `Result_aiy_aiz_negative.json` の遺伝子（インデックス 0）の解析データ。AIZ–SMB シナプスを段階的に弱めたもの。
* **ファイル**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` ～ `n_vs_c_10.txt`

#### 4. `Result_aiz_smb/`

* **内容**: `concentration_memory/Result_aiy_aiz_negative_0.json` の遺伝子（インデックス 0, 9, 15）の解析データ。AIZ–SMB および SMB–SMB シナプスを 0.9 倍にスケーリング。
* **ファイル**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` ～ `n_vs_c_2.txt`
* **遺伝子インデックスと条件**:

  * `0`: 高塩条件で飼育
  * `1`: 中塩条件
  * `2`: 低塩条件

---

### 🔸 `bias/`

#### 5. `Result_smb_0/`

* **内容**: `Result_aiy_aiz_negative.json` の遺伝子（インデックス 0）の解析データ。SMB ニューロンのバイアスを段階的に増加。
* **ファイル**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` ～ `n_vs_c_10.txt`

#### 6. `Result_smb/`

* **内容**: `concentration_memory/Result_aiy_aiz_negative_0.json` の遺伝子（インデックス 0, 9, 15）の解析データ。SMB バイアスを一律 –0.05 に設定。
* **ファイル**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` ～ `n_vs_c_2.txt`
* **備考**: 遺伝子インデックスと条件は `synapse` サブディレクトリと同じ。

---

## 🔚 備考

* **全ファイルのカラム形式**:

  1. 入力変数（例：方位、勾配）
  2. 屈曲率
  3. 統計指標（標準偏差、最大値、最小値など）

* すべてプレーンテキスト形式であり、Python や Rust などのスクリプトで容易に読み込み可能。

* ディレクトリ構造により、どの遺伝子・条件で得られたデータかが一目で分かる。

* `gene/` ディレクトリと併用することで、シナプス操作が線虫の行動に与える影響を体系的に評価可能。
