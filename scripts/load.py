import toml
import json
import re


def load_simulation_setting_toml(file_path):
    """
    file_path: 読み込むシミュレーション設定に関するtomlファイルのパス
    """
    with open(file_path, "r") as toml_file:
        data = toml.load(toml_file)
    return data


def load_result_json(file_path):
    """
    file_path: 読み込む遺伝子に関するjsonファイルのパス
    """
    with open(file_path, "r") as json_file:
        data = json.load(json_file)
    return data


def load_output_txt(file_path):
    """
    file_path: 読み込むアウトプットに関するtxtファイルのパス

    データは列で読み込む（カンマおよびスペースで分割）
    """
    with open(file_path, "r") as txt_file:
        # ファイルの各行を読み込んでリストに格納
        lines = txt_file.read().split("\n")

    # 「#」で始まる行を除外
    lines = [line for line in lines if not line.startswith("#")]

    # 分割パターンの正規表現
    split_pattern = re.compile(r"[,\s]+")

    # 各行での最大列数を求める
    max_columns = max(
        len([value for value in split_pattern.split(line) if value.strip()])
        for line in lines
        if line
    )

    # 列ごとのデータを格納するためのリストを初期化
    columns = [[] for _ in range(max_columns)]

    # 行ごとにデータを列に振り分ける
    for line in lines:
        if line:
            # 各行のデータを正規表現で分割し、空でない文字列を除外してからfloatに変換
            values = [
                float(value) for value in split_pattern.split(line) if value.strip()
            ]
            # 各列にデータを振り分け
            for i in range(len(values)):
                columns[i].append(values[i])

    return columns
