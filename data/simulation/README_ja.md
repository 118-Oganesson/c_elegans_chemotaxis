# Description of the `simulation` Directory

This directory contains multiple text files with the results of behavioral simulations of *C. elegans* based on genetic parameter sets. Each file records statistical analyses of various behavioral metrics (e.g., curving rate) computed under different conditions. Files are organized by name and subdirectory.

---

## 📁 File & Directory Overview

### 1. `with_constraining_AIY-AIZ/`

* **Contents**: Simulation results for the best-performing individual whose genes were optimized under an inhibitory constraint on the AIY–AIZ synapse.
* **Files**:

  * `bearing_vs_curving_rate_0.txt`
  * `normal_gradient_vs_curving_rate_0.txt`
  * `translational_gradient_vs_curving_rate_0.txt`
* **Columns**:

  * **bearing\_vs\_curving\_rate**: `[bearing, curving rate, std, max, min]`
  * **normal\_gradient\_vs\_curving\_rate**: `[normal gradient, curving rate, std, max, min]`
  * **translational\_gradient\_vs\_curving\_rate**: `[translational gradient, curving rate, std, turning bias (+), std (+), turning bias (–), std (–)]`

---

## 📁 `concentration_memory/`

Holds analysis results for chemotaxis simulations that incorporate salt-concentration memory.

### 2. `Result_aiy_aiz_negative_0/`

* **Contents**: Analysis data for the top gene (index 0) from `Result_aiy_aiz_negative.json`, with stepwise modification of the ASER–AIY synapse.
* **Files**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` through `n_vs_c_15.txt`

---

## 📁 `starvation/`

Contains simulation analyses under starvation conditions.

### 🔸 `synapse/`

#### 3. `Result_aiz_smb_0/`

* **Contents**: Analysis data for the 0th gene in `Result_aiy_aiz_negative.json`, with stepwise weakening of the AIZ–SMB synapse.
* **Files**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` through `n_vs_c_10.txt`

#### 4. `Result_aiz_smb/`

* **Contents**: Analysis data for gene indices 0, 9, and 15 from `concentration_memory/Result_aiy_aiz_negative_0.json`, scaling both AIZ–SMB and SMB–SMB synapses by 0.9×.
* **Files**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` through `n_vs_c_2.txt`
* **Notes on Gene Indices**:

  * `0`: High-salt rearing
  * `1`: Medium
  * `2`: Low-salt rearing

### 🔸 `bias/`

#### 5. `Result_smb_0/`

* **Contents**: Analysis data for the 0th gene in `Result_aiy_aiz_negative.json`, with stepwise increases in SMB neuron bias.
* **Files**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` through `n_vs_c_10.txt`

#### 6. `Result_smb/`

* **Contents**: Analysis data for gene indices 0, 9, and 15 from `concentration_memory/Result_aiy_aiz_negative_0.json`, setting SMB bias uniformly to –0.05.
* **Files**:

  * `normal_gradient_vs_curving_rate/n_vs_c_0.txt` through `n_vs_c_2.txt`
* **Notes**: Gene indices and conditions match those in the `synapse` subdirectory.

---

## 🔚 Notes

* **Column format (all files)**:

  1. Input variable (e.g., bearing, gradient)
  2. Curving rate
  3. Statistical metrics (std, max, min, etc.)

* All files are plain text and can be readily loaded by scripts (Python, Rust, etc.).

* Directory structure makes it easy to see at a glance which gene and condition produced each data set.

* When used alongside the `gene/` directory, these results allow systematic evaluation of how synaptic modulations affect worm behavior.
