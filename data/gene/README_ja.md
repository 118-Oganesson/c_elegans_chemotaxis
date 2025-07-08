# gene ディレクトリの説明

このディレクトリには、遺伝的アルゴリズムによって最適化された遺伝子データとその評価値を含む複数の JSON ファイルが格納されています。各ファイルは異なる実験条件下で最適化されたデータを記録しており、ファイル名やサブディレクトリによって分類されています。

---

## 📁 ファイル・ディレクトリ一覧と説明

### 1. `Result.json`

- **内容**：特別な制約なしで最適化された遺伝子データ。
- **構成**：
  - `value`: 評価値
  - `gene`: 遺伝子の数値リスト（シナプス強度など）

```json
[
    {
        "value": 0.0,
        "gene": [
            -0.8094022576319283,
            -0.6771492613425638,
            ...
        ]
    },
    {
        "value": 0.1,
        "gene": [
            -0.8094022576319283,
            -0.6771492613425638,
            ...
        ]
    }
]
```

---

### 2. `Result_aiy_aiz_negative.json`

- **内容**：AIY-AIZ のシナプスを抑制的に固定した状態で最適化した遺伝子データ。
- **構成**：`Result.json` と同様

---

## 📁 concentration_memory/

条件付きの記憶保持機構に関する最適化や変化を扱うデータ群です。

### 3. `Result_aiy_aiz_negative_0.json`

- **内容**：`Result_aiy_aiz_negative.json` の 0 番目の遺伝子を基に、ASER-AIY シナプスを徐々に興奮性に変化させたデータ。
- **value**：ASER-AIY のシナプス特性に加える値（大きいほど興奮性）

### 4. `Result_aiy_aiz_negative_1.json`

- **内容**：`Result_aiy_aiz_negative.json` の 1 番目の遺伝子を基に、同様の条件でシナプス特性の変化を調べたデータ。
- **value**：ASER-AIY のシナプス特性に加える値（大きいほど興奮性）

---

## 📁 starvation/

飢餓条件下における変化を反映した遺伝子データ。

### 🔸 synapse/

### 5. `Result_aiz_smb_0.json`

- **内容**：`Result_aiy_aiz_negative.json` の 0 番目の遺伝子を基に、AIZ-SMB のシナプスを徐々に弱めたデータ。
- **value**：シナプスにかけるスケーリング係数

### 6. `Result_aiz_smb.json`

- **内容**：`concentration_memory/Result_aiy_aiz_negative_0.json` の 0, 9, 15 番目の遺伝子を基に、AIZ-SMB および SMB-SMB のシナプスを 0.9 倍にスケーリング。
- **備考**：遺伝子番号
  - `0`: 高塩濃度育成
  - `1`: 中間
  - `2`: 低塩濃度育成

---

### 🔸 bias/

### 7. `Result_smb_0.json`

- **内容**：`Result_aiy_aiz_negative.json` の 0 番目の遺伝子を基に、SMB のバイアスの影響を徐々に強めたデータ。
- **value**：SMBのバイアスから引く値（大きいほど影響大）

### 8. `Result_smb.json`

- **内容**：`concentration_memory/Result_aiy_aiz_negative_0.json` の 0, 9, 15 番目の遺伝子に対して、SMB バイアスを一律 -0.05 したデータ。
- **備考**：遺伝子番号
  - `0`: 高塩濃度育成
  - `1`: 中間
  - `2`: 低塩濃度育成

---

## 🔚 備考

- 全ファイルに共通して、`gene` フィールドは最適化された遺伝子パラメータのリストです。
- `value` フィールドは、最適化中に用いられた評価関数の値、または操作パラメータを示します。
- JSON データは Python や Rust などのスクリプトから簡単に読み込める形式です。
