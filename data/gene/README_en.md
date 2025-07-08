Here is a natural English translation of your explanation for the **gene** directory:

---

# Description of the `gene` Directory

This directory contains multiple JSON files storing genetic data optimized by a genetic algorithm, along with their evaluation values. Each file corresponds to optimization results under different experimental conditions, organized by file names and subdirectories.

---

## 📁 Files and Directories Overview

### 1. `Result.json`

* **Contents:** Genetic data optimized without special constraints.
* **Structure:**

  * `value`: Evaluation score
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

* **Contents:** Genetic data optimized with the AIY-AIZ synapse fixed as inhibitory.
* **Structure:** Same as `Result.json`

---

## 📁 `concentration_memory/`

Data related to optimization and changes involving conditional memory retention mechanisms.

### 3. `Result_aiy_aiz_negative_0.json`

* **Contents:** Based on the 0th gene from `Result_aiy_aiz_negative.json`, this dataset gradually changes the ASER-AIY synapse to be more excitatory.
* **`value`:** The parameter added to the ASER-AIY synapse property (larger values mean more excitatory)

### 4. `Result_aiy_aiz_negative_1.json`

* **Contents:** Similar to the above, but based on the 1st gene from `Result_aiy_aiz_negative.json`.
* **`value`:** The parameter added to the ASER-AIY synapse property (larger values mean more excitatory)

---

## 📁 `starvation/`

Genetic data reflecting changes under starvation conditions.

### 🔸 `synapse/`

### 5. `Result_aiz_smb_0.json`

* **Contents:** Based on the 0th gene from `Result_aiy_aiz_negative.json`, this data gradually weakens the AIZ-SMB synapse.
* **`value`:** Scaling factor applied to the synapse strength

### 6. `Result_aiz_smb.json`

* **Contents:** Based on the 0th, 9th, and 15th genes from `concentration_memory/Result_aiy_aiz_negative_0.json`, synapses AIZ-SMB and SMB-SMB are scaled by 0.9.
* **Notes:** Gene indices correspond to different salt concentration conditions during development:

  * `0`: High salt concentration
  * `1`: Medium
  * `2`: Low salt concentration

---

### 🔸 `bias/`

### 7. `Result_smb_0.json`

* **Contents:** Based on the 0th gene from `Result_aiy_aiz_negative.json`, this dataset gradually increases the influence of SMB bias.
* **`value`:** Amount subtracted from the SMB bias (larger values mean stronger effect)

### 8. `Result_smb.json`

* **Contents:** For the 0th, 9th, and 15th genes from `concentration_memory/Result_aiy_aiz_negative_0.json`, SMB bias is uniformly decreased by -0.05.
* **Notes:** Gene indices correspond to different salt concentration conditions during development:

  * `0`: High salt concentration
  * `1`: Medium
  * `2`: Low salt concentration

---

## 🔚 Notes

* In all files, the `gene` field contains the list of optimized gene parameters.
* The `value` field represents either the evaluation function value during optimization or a parameter controlling modifications.
* The JSON format can be easily loaded by Python, Rust, or other scripts for analysis.
