# Description of the `gene` Directory

This directory contains multiple JSON files, each holding genetic parameter data optimized by a genetic algorithm along with their associated fitness values. Each file records data optimized under different experimental conditions, and is organized by filename or subdirectory.

---

## 📁 File & Directory Overview

### 1. `Result.json`

* **Contents**: Genetic data optimized with no special constraints.
* **Structure**:

  * `value`: Fitness score
  * `gene`: List of numerical gene parameters (e.g., synaptic weights)

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

* **Contents**: Genetic data optimized under an inhibitory constraint on the AIY–AIZ synapse.
* **Structure**: Same as `Result.json`.

---

## 📁 `concentration_memory/`

Contains data relating to optimization and manipulation of chemotaxis based on salt-concentration memory.

### 3. `Result_aiy_aiz_negative_0.json`

* **Contents**: Starting from the 0th entry in `Result_aiy_aiz_negative.json`, this file gradually shifts the ASER–AIY synapse toward excitatory behavior.
* **`value`**: Amount added to the ASER–AIY synapse parameter (larger values = stronger excitation).

### 4. `Result_aiy_aiz_negative_1.json`

* **Contents**: Similar to #3, but based on the 1st entry of `Result_aiy_aiz_negative.json`.
* **`value`**: Amount added to the ASER–AIY synapse parameter.

---

## 📁 `starvation/`

Holds genetic data reflecting changes under starvation conditions.

### 🔸 `synapse/`

#### 5. `Result_aiz_smb_0.json`

* **Contents**: Based on the 0th entry of `Result_aiy_aiz_negative.json`, this file gradually weakens the AIZ–SMB synapse.
* **`value`**: Scaling factor applied to that synapse.

#### 6. `Result_aiz_smb.json`

* **Contents**: From `concentration_memory/Result_aiy_aiz_negative_0.json`, uses entries 0, 9, and 15 to scale AIZ–SMB and SMB–SMB synapses by 0.9×.
* **Notes**:

  * Gene index `0`: High-salt cultivation
  * Gene index `1`: Medium
  * Gene index `2`: Low-salt cultivation

### 🔸 `bias/`

#### 7. `Result_smb_0.json`

* **Contents**: Based on the 0th entry of `Result_aiy_aiz_negative.json`, this file incrementally increases the influence of the SMB neuron’s bias.
* **`value`**: Amount subtracted from the SMB bias (larger values = greater effect).

#### 8. `Result_smb.json`

* **Contents**: From `concentration_memory/Result_aiy_aiz_negative_0.json`, uses entries 0, 9, and 15, setting the SMB bias uniformly to –0.05.
* **Notes**: Gene indices are the same as in `Result_aiz_smb.json`.

---

## 🔚 Notes

* In every file, the `gene` field is the list of optimized genetic parameters.
* The `value` field represents either the fitness score from the optimization or the parameter used in the manipulation.
* All JSON files can be easily loaded by Python, Rust, or similar scripts.
