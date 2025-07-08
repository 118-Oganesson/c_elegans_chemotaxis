from scripts import load
from scripts import oed
import os
import json
import numpy as np


def gene_aser_aiy_synapse_change(
    gene_number, start, stop, num, gene_path="../result/Result_aiy_aiz_negative.json"
):
    """
    ASER-AIYのシナプス特性を変化させる。
    gene_number: 改変する遺伝子の番号
    start: 10,11の結合(ASER-AIY)に対して加える範囲の最小
    stop: 10,11の結合(ASER-AIY)に対して加える範囲の最大
    num: 範囲の分割数
    gene_path: 遺伝子が保存されているJSONファイルのpath
    """

    result = load.load_result_json(gene_path)
    base_gene = np.array(result[gene_number]["gene"])

    change_result = []
    for delta_gene in np.linspace(start, stop, num + 1):
        zero_gene = np.zeros_like(base_gene)
        zero_gene[10:12] = [delta_gene] * 2
        change_gene = base_gene + zero_gene
        change_gene = np.clip(change_gene, -1, 1)
        change_result.append({"value": delta_gene, "gene": change_gene.tolist()})

    json_data = json.dumps(change_result, indent=1)
    os.makedirs("../result/concentration_memory", exist_ok=True)
    with open(
        "../result/concentration_memory/Result_aiy_aiz_negative_{}.json".format(
            gene_number
        ),
        "w",
    ) as json_file:
        json_file.write(json_data)


def gene_ase_aiy_synapse_change(
    gene_number, start, stop, num, gene_path="../result/Result_aiy_aiz_negative.json"
):
    """
    ASER-AIYのシナプス特性を変化させる。
    gene_number: 改変する遺伝子の番号
    start: 8,9,10,11の結合(ASE-AIY)に対して加える範囲の最小
    stop: 8,9,10,11の結合(ASE-AIY)に対して加える範囲の最大
    num: 範囲の分割数
    gene_path: 遺伝子が保存されているJSONファイルのpath
    """

    result = load.load_result_json(gene_path)
    base_gene = np.array(result[gene_number]["gene"])

    change_result = []
    for delta_gene in np.linspace(start, stop, num + 1):
        zero_gene = np.zeros_like(base_gene)
        zero_gene[8:12] = [delta_gene] * 4
        change_gene = base_gene + zero_gene
        change_gene = np.clip(change_gene, -1, 1)
        change_result.append({"value": delta_gene, "gene": change_gene.tolist()})

    json_data = json.dumps(change_result, indent=1)
    os.makedirs("../result/concentration_memory", exist_ok=True)
    with open(
        "../result/concentration_memory/Result_aiy_aiz_negative_{}_ase_aiy_positive.json".format(
            gene_number
        ),
        "w",
    ) as json_file:
        json_file.write(json_data)


def gene_asel_aiy_delete(
    gene_number,
    gene_path=[
        "../result/Result.json",
        "../result/Result_aiy_aiz_negative.json",
        "../result/concentration_memory/Result_aiy_aiz_negative_0_optimize_aser_aiy_positive.json",
    ],
):
    """
    ASELの機能を停止させる。(8,9の結合を切る)
    gene_number: 改変する遺伝子の番号
    gene_path: 遺伝子が保存されているJSONファイルのpathのリスト[制約無しモデル，制約ありモデル，低塩濃度モデル]
    """

    base_gene = []
    for path in gene_path:
        result = load.load_result_json(path)
        base_gene.append(result[gene_number]["gene"])

    change_result = []
    for gene in base_gene:
        gene[8:10] = [0.0] * 2
        change_result.append({"value": 0.0, "gene": gene})

    json_data = json.dumps(change_result, indent=1)
    os.makedirs("../result/newron_deletion", exist_ok=True)
    with open(
        "../result/newron_deletion/Result_asel_{}.json".format(gene_number),
        "w",
    ) as json_file:
        json_file.write(json_data)


def gene_aiz_smb_synapse_change(
    gene_number=[0, 9, 15],
    scaling=0.9,
    gene_path="../result/concentration_memory/Result_aiy_aiz_negative_0.json",
):
    """
    AIZ-SMBのシナプス特性を変化させる。
    gene_number: 改変する遺伝子の番号[高塩濃度, 中間, 低塩濃度]
    scaling: 14,15の結合(AIZ-SMB)に対してかける値
    gene_path: 遺伝子が保存されているJSONファイルのpath
    """

    result = load.load_result_json(gene_path)
    base_gene = []
    for i in gene_number:
        base_gene.append(result[i]["gene"])

    change_result = []
    for gene in base_gene:
        gene[14] *= scaling
        gene[15] *= scaling
        gene[16] *= scaling
        gene[17] *= scaling
        change_result.append({"value": 0.0, "gene": gene})

    json_data = json.dumps(change_result, indent=1)
    os.makedirs("../result/starvation/synapse", exist_ok=True)
    with open(
        "../result/starvation/synapse/Result_aiz_smb.json",
        "w",
    ) as json_file:
        json_file.write(json_data)


def gene_aiz_smb_synapse_gradually_change(
    gene_number, start, stop, num, gene_path="../result/Result_aiy_aiz_negative.json"
):
    """
    AIZ-SMBのシナプス特性を変化させる。
    gene_number: 改変する遺伝子の番号
    start: 14,15の結合(AIZ-SMB)に対してかける範囲の最小
    stop: 14,15の結合(AIZ-SMB)に対してかける範囲の最大
    num: 範囲の分割数
    gene_path: 遺伝子が保存されているJSONファイルのpath
    """

    result = load.load_result_json(gene_path)
    base_gene = np.array(result[gene_number]["gene"])

    change_result = []
    for delta_gene in np.linspace(start, stop, num + 1):
        one_gene = np.ones_like(base_gene)
        one_gene[14:16] = [delta_gene] * 2
        change_gene = base_gene * one_gene
        change_gene = np.clip(change_gene, -1, 1)
        change_result.append({"value": delta_gene, "gene": change_gene.tolist()})

    json_data = json.dumps(change_result, indent=1)
    os.makedirs("../result/starvation/synapse", exist_ok=True)
    with open(
        "../result/starvation/synapse/Result_aiz_smb_negative_{}.json".format(
            gene_number
        ),
        "w",
    ) as json_file:
        json_file.write(json_data)


def gene_smb_bias_change(
    gene_number=[0, 9, 15],
    scaling=0.05,
    gene_path="../result/concentration_memory/Result_aiy_aiz_negative_0.json",
):
    """
    SMBのバイアスを変化させる。
    gene_number: 改変する遺伝子の番号[高塩濃度, 中間, 低塩濃度]
    scaling: 6,7のバイアス(SMB)に対して引く値
    gene_path: 遺伝子が保存されているJSONファイルのpath
    """
    result = load.load_result_json(gene_path)
    base_gene = []
    for i in gene_number:
        base_gene.append(result[i]["gene"])

    change_result = []
    for gene in base_gene:
        gene[6] -= scaling
        gene[7] -= scaling
        change_result.append({"value": 0.0, "gene": gene})

    json_data = json.dumps(change_result, indent=1)
    os.makedirs("../result/starvation/bias", exist_ok=True)
    with open(
        "../result/starvation/bias/Result_smb.json",
        "w",
    ) as json_file:
        json_file.write(json_data)


def gene_smb_bias_gradually_change(
    gene_number, start, stop, num, gene_path="../result/Result_aiy_aiz_negative.json"
):
    """
    SMBのバイアスを変化させる。
    gene_number: 改変する遺伝子の番号
    start: 6,7のバイアス(SMB)に対して引く範囲の最小
    stop: 6,7のバイアス(SMB)に対して引く範囲の最大
    num: 範囲の分割数
    gene_path: 遺伝子が保存されているJSONファイルのpath
    """

    result = load.load_result_json(gene_path)
    base_gene = np.array(result[gene_number]["gene"])

    change_result = []
    for delta_gene in np.linspace(start, stop, num + 1):
        zero_gene = np.zeros_like(base_gene)
        zero_gene[6:8] = [delta_gene] * 2
        change_gene = base_gene - zero_gene
        change_gene = np.clip(change_gene, -1, 1)
        change_result.append({"value": delta_gene, "gene": change_gene.tolist()})

    json_data = json.dumps(change_result, indent=1)
    os.makedirs("../result/starvation/bias", exist_ok=True)
    with open(
        "../result/starvation/bias/Result_smb_{}.json".format(gene_number),
        "w",
    ) as json_file:
        json_file.write(json_data)


def parameter_output_latex(gene, decimal_places):
    N, M, theta, w_on, w_off, w, g, w_osc, w_nmj = oed.weight(gene)
    N, M, theta, w_on, w_off, w, g, w_osc, w_nmj = [
        np.round(var, decimal_places)
        for var in [N, M, theta, w_on, w_off, w, g, w_osc, w_nmj]
    ]

    name = [
        "AIYL",
        "AIYR",
        "AIZL",
        "AIZR",
        "SMBVL",
        "SMBDL",
        "SMBDR",
        "SMBVR",
    ]

    latex_output = []

    latex_output.append("\\begin{align*}")
    latex_output.append(f"N &= {N:.{decimal_places}f} \\\\")
    latex_output.append(f"M &= {M:.{decimal_places}f} \\\\")

    for i, data in enumerate(theta):
        latex_output.append(
            f"\\theta_{{\\mathrm{{{name[i]}}}}} &= {data:.{decimal_places}f} \\\\"
        )

    for i, data in enumerate(w_on):
        if data != 0.0:
            latex_output.append(
                f"w_{{\\mathrm{{{name[i]}}}}}^{{\\mathrm{{ON}}}} &= {data:.{decimal_places}f} \\\\"
            )

    for i, data in enumerate(w_off):
        if data != 0.0:
            latex_output.append(
                f"w_{{\\mathrm{{{name[i]}}}}}^{{\\mathrm{{OFF}}}} &= {data:.{decimal_places}f} \\\\"
            )

    for i in range(8):
        for j in range(8):
            if w[i][j] != 0.0 and i != j:
                latex_output.append(
                    f"w_{{\\mathrm{{{name[i]}}}, \\mathrm{{{name[j]}}}}} &= {w[i][j]:.{decimal_places}f} \\\\"
                )

    for i in range(8):
        for j in range(8):
            if w[i][j] != 0.0 and i == j:
                latex_output.append(
                    f"w_{{\\mathrm{{{name[i]}}}}}^{{\\mathrm{{self}}}} &= {w[i][j]:.{decimal_places}f} \\\\"
                )

    for i in range(8):
        for j in range(8):
            if g[i][j] != 0.0 and i != j:
                latex_output.append(
                    f"g_{{\\mathrm{{{name[i]}}}, \\mathrm{{{name[j]}}}}} &= {g[i][j]:.{decimal_places}f} \\\\"
                )

    latex_output.append(f"w_{{\\mathrm{{OSC}}}} &= {w_osc[4]:.{decimal_places}f} \\\\")
    latex_output.append(f"w_{{\\mathrm{{NMJ}}}} &= {w_nmj:.{decimal_places}f} \\\\")
    latex_output.append("\\end{align*}")

    return "\n".join(latex_output)
