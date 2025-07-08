import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import animation
from PIL import Image, ImageOps
from scripts import oed
from scripts import figure
import time as tm


def calculate_trajectory(gene_angle_list):
    gene, angle, c_mode = gene_angle_list
    return oed.klinotaxis_animation(gene, angle, c_mode)


def compute_concentration_map(x_range, y_range, c_mode):
    """指定されたxおよびy範囲に対する濃度関数を計算します。"""
    alpha, x_peak, y_peak, dt, T, f, v, time, tau, c_0, lambda_ = oed.constant(
        "setting_animation"
    )
    x = np.linspace(x_range[0], x_range[1], num=500)
    y = np.linspace(y_range[0], y_range[1], num=500)
    x, y = np.meshgrid(x, y)
    if c_mode == 0:
        z = oed.c_(alpha, x, y, x_peak, y_peak)
    elif c_mode == 1:
        z = oed.c_gauss(c_0, lambda_, x, y, x_peak, y_peak)
    elif c_mode == 2:
        z = oed.c_two_gauss(c_0, lambda_, x, y, x_peak, y_peak)
    return z


def create_light_colormap():
    """明るく控えめなカラーマップを作成します。"""
    colors = [
        "#ffffff",
        "#f0f8ff",
        "#e0f0ff",
        "#c0e0ff",
        "#a0d0ff",
        "#80c0ff",
        "#60b0ff",
        "#40a0ff",
        "#2090ff",
    ]
    cmap = LinearSegmentedColormap.from_list("natural_blue_gradient", colors, N=256)
    return cmap


def single_trajectory_animation(
    gene,
    c_mode=1,
    lines_number=10,
    out_file_path="single_output.mp4",
    skip=10,
    downsample_factor=10,
    smooth_factor=0.1,
    padding=1,
    figsize=(10, 4),
    fontsize=14,
):
    """
    遺伝子の軌跡アニメーションを作成し、保存します。
    使用例：
    result = load.load_result_json("../result/Result_aiy_aiz_negative.json")
    gene = result[0]["gene"]
    animation.single_trajectory_animation(gene)
    """
    alpha, x_peak, y_peak, dt, T, f, v, time, tau, c_0, lambda_ = oed.constant(
        "setting_animation"
    )
    fig, ax = plt.subplots(figsize=figsize)

    # 並列計算のためのプロセス数を決定
    process = min(lines_number, multiprocessing.cpu_count())
    gene_angle_list = [
        [gene, angle, c_mode]
        for angle in np.arange(0, 2 * np.pi, 2 * np.pi / lines_number)
    ]

    # 遺伝子の軌跡を並列で計算
    with multiprocessing.Pool(process) as pool:
        results = pool.map(calculate_trajectory, gene_angle_list)

    # 結果をダウンサンプル
    results = [[b[::downsample_factor] for b in a] for a in results]
    all_lines = []
    for r in results:
        lines = figure.single_line_stacks(r[0], r[1])
        all_lines.append(lines)

    # 広い範囲の濃度マップを計算
    broad_x_range = [-15, 15]
    broad_y_range = [-15, 15]
    broad_concentration_map = compute_concentration_map(
        broad_x_range, broad_y_range, c_mode
    )

    # 背景画像をロードして準備
    img = Image.open("./c_elegans.png").convert("RGBA")
    img_rotated = img.rotate(-60, expand=True)
    img_resized = img_rotated.resize((300, 300))
    img_flipped = ImageOps.mirror(img_resized)
    img_resized = np.array(img_resized)
    img_flipped = np.array(img_flipped)

    def get_imagebox(forward):
        return OffsetImage(img_resized if forward else img_flipped, zoom=0.1)

    annotation_boxes = [
        AnnotationBbox(
            get_imagebox(True),
            (0, 0),
            xybox=(0, 0),
            xycoords="data",
            boxcoords="offset points",
            frameon=False,
        )
        for _ in range(lines_number)
    ]
    for ab in annotation_boxes:
        ax.add_artist(ab)

    # 開始点とピーク点を定義
    starting_point = [0, 0]
    peak = [x_peak, y_peak]

    # 開始点とピーク点をプロット
    ax.scatter(*starting_point, s=15, color="black")
    ax.scatter(*peak, s=15, color="black")
    y_max = 1
    ax.vlines(
        starting_point[0],
        starting_point[1],
        y_max,
        color="black",
        linestyle="-",
        linewidth=0.5,
    )
    ax.vlines(peak[0], peak[1], y_max, color="black", linestyle="-", linewidth=0.5)
    ax.text(
        starting_point[0],
        y_max + 0.1,
        "Starting Point",
        horizontalalignment="center",
        fontsize=fontsize,
    )
    ax.text(
        peak[0],
        y_max + 0.1,
        "Gradient Peak",
        horizontalalignment="center",
        fontsize=fontsize,
    )

    ax.axis("off")
    ax.set_aspect("equal")
    ax.text(4.5, -0.95, "1 cm", horizontalalignment="center", fontsize=fontsize)
    ax.hlines(-1, 4, 5, color="black", linestyle="-", linewidth=1.5)

    # アニメーション用のラインコレクションを初期化
    line_collections = [
        LineCollection([], colors="gray", linewidth=1) for _ in range(lines_number)
    ]
    for lc in line_collections:
        ax.add_collection(lc)

    prev_lengths = [0] * lines_number
    prev_min_x, prev_max_x, prev_min_y, prev_max_y = None, None, None, None

    # 経過時間を表示するテキスト
    time_text = fig.text(
        0.4,
        0.85,
        "",
        horizontalalignment="center",
        fontsize=fontsize * 1.5,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),
    )

    def animate(frame):
        """プロットを更新するためのアニメーション関数。"""
        nonlocal prev_min_x, prev_max_x, prev_min_y, prev_max_y

        def update_lines(
            all_lines,
            line_collections,
            annotation_boxes,
            prev_lengths,
            ax,
            starting_point,
            peak,
            prev_min_x,
            prev_max_x,
            prev_min_y,
            prev_max_y,
        ):
            """ラインを更新し、軸の範囲を調整します。"""
            min_x = min(starting_point[0], peak[0])
            max_x = max(starting_point[0], peak[0])
            min_y = min(starting_point[1], peak[1])
            max_y = max(starting_point[1], peak[1])

            for i, lines in enumerate(all_lines):
                prev_length = prev_lengths[i]
                current_length = min(len(lines), frame * skip + 1)
                new_segments = lines[prev_length:current_length]
                if new_segments:
                    segments = list(line_collections[i].get_segments())
                    segments.extend(new_segments)
                    line_collections[i].set_segments(segments)

                if current_length > 0:
                    latest_point = lines[current_length - 1][-1]
                    annotation_boxes[i].xy = (latest_point[0], latest_point[1])

                    if current_length > 1:
                        previous_point = lines[current_length - 2][-1]
                        forward = latest_point[0] >= previous_point[0]
                        annotation_boxes[i].offsetbox = get_imagebox(forward)

                    min_x = min(min_x, latest_point[0])
                    max_x = max(max_x, latest_point[0])
                    min_y = min(min_y, latest_point[1])
                    max_y = max(max_y, latest_point[1])

                prev_lengths[i] = current_length

            target_min_x = min_x - padding
            target_max_x = max_x + padding
            target_min_y = min_y - padding
            target_max_y = max_y + padding

            if prev_min_x is None:
                prev_min_x, prev_max_x, prev_min_y, prev_max_y = (
                    target_min_x,
                    target_max_x,
                    target_min_y,
                    target_max_y,
                )
            else:
                prev_min_x += smooth_factor * (target_min_x - prev_min_x)
                prev_max_x += smooth_factor * (target_max_x - prev_max_x)
                prev_min_y += smooth_factor * (target_min_y - prev_min_y)
                prev_max_y += smooth_factor * (target_max_y - prev_max_y)

            ax.set_xlim(prev_min_x, prev_max_x)
            ax.set_ylim(prev_min_y, prev_max_y)

            return prev_min_x, prev_max_x, prev_min_y, prev_max_y

        prev_min_x, prev_max_x, prev_min_y, prev_max_y = update_lines(
            all_lines,
            line_collections,
            annotation_boxes,
            prev_lengths,
            ax,
            starting_point,
            peak,
            prev_min_x,
            prev_max_x,
            prev_min_y,
            prev_max_y,
        )

        # 必要に応じて濃度マップを更新
        if (
            prev_min_x < broad_x_range[0]
            or prev_max_x > broad_x_range[1]
            or prev_min_y < broad_y_range[0]
            or prev_max_y > broad_y_range[1]
        ):
            z = compute_concentration_map(
                [prev_min_x, prev_max_x], [prev_min_y, prev_max_y], c_mode
            )
            im.set_data(z)

        # 経過時間を更新
        elapsed_time = frame * (time / (max(len(lines) for lines in all_lines) // skip))
        time_text.set_text(f"Time: {int(elapsed_time)} s")

    # 総フレーム数を計算
    total_frames = max(len(lines) for lines in all_lines) // skip

    # 濃度マップとカラーバーの作成
    cmap = create_light_colormap()
    im = ax.imshow(
        broad_concentration_map,
        cmap=cmap,
        origin="lower",
        extent=[broad_x_range[0], broad_x_range[1], broad_y_range[0], broad_y_range[1]],
        alpha=0.5,
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label("Concentration (mM)", fontsize=fontsize)

    # アニメーションの作成
    ani = animation.FuncAnimation(
        fig, animate, frames=total_frames, interval=100, repeat=True
    )

    # アニメーションを保存
    start_time = tm.time()
    ani.save(out_file_path, writer="ffmpeg", dpi=300)
    save_time = tm.time() - start_time
    print(f"Time taken to save the animation: {save_time:.2f} seconds")


def dual_trajectory_animation(
    gene1,
    gene2,
    c_mode=1,
    lines_number=10,
    out_file_path="dual_output.mp4",
    skip=10,
    downsample_factor=10,
    smooth_factor=0.1,
    padding=1,
    figsize=(10, 8),
    fontsize=14,
):
    """
    遺伝子の軌跡アニメーションを作成し、保存します。
    使用例：
    result_1 = load.load_result_json("../result/Result_aiy_aiz_negative.json")
    result_2 = load.load_result_json(
        "../result/concentration_memory/Result_aiy_aiz_negative_0.json"
    )
    gene_1 = result_1[0]["gene"]
    gene_2 = result_2[15]["gene"]
    animation.dual_trajectory_animation(gene_1, gene_2)
    """
    # 定数の取得
    alpha, x_peak, y_peak, dt, T, f, v, time, tau, c_0, lambda_ = oed.constant(
        "setting_animation"
    )

    # マルチプロセスのためのプロセス数設定
    process = min(lines_number, multiprocessing.cpu_count())

    def calculate_lines(gene):
        # 遺伝子の角度ごとのトラジェクトリを計算
        gene_angle_list = [
            [gene, angle, c_mode]
            for angle in np.arange(0, 2 * np.pi, 2 * np.pi / lines_number)
        ]

        with multiprocessing.Pool(process) as pool:
            results = pool.map(calculate_trajectory, gene_angle_list)

        # データの間引き
        results = [[b[::downsample_factor] for b in a] for a in results]

        all_lines = [figure.single_line_stacks(r[0], r[1]) for r in results]
        return all_lines

    all_lines1 = calculate_lines(gene1)
    all_lines2 = calculate_lines(gene2)

    # 広い範囲の濃度マップを計算
    broad_x_range = [-15, 15]
    broad_y_range = [-15, 15]
    broad_concentration_map = compute_concentration_map(
        broad_x_range, broad_y_range, c_mode
    )

    # C. elegansの画像読み込みと加工
    img = Image.open("./c_elegans.png").convert("RGBA")
    img_rotated = img.rotate(-60, expand=True)
    img_resized = img_rotated.resize((300, 300))
    img_flipped = ImageOps.mirror(img_resized)
    img_resized = np.array(img_resized)
    img_flipped = np.array(img_flipped)

    def get_imagebox(forward):
        """アニメーションの進行方向に応じて画像を取得"""
        return OffsetImage(img_resized if forward else img_flipped, zoom=0.1)

    # グラフの作成
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

    def setup_axis(ax):
        """各軸の設定"""
        annotation_boxes = [
            AnnotationBbox(
                get_imagebox(True),
                (0, 0),
                xybox=(0, 0),
                xycoords="data",
                boxcoords="offset points",
                frameon=False,
            )
            for _ in range(lines_number)
        ]
        for ab in annotation_boxes:
            ax.add_artist(ab)

        starting_point = [0, 0]
        peak = [x_peak, y_peak]

        # 開始点とピークの設定
        ax.scatter(*starting_point, s=15, color="black")
        ax.scatter(*peak, s=15, color="black")
        ax.vlines(
            starting_point[0],
            starting_point[1],
            1,
            color="black",
            linestyle="-",
            linewidth=0.5,
        )
        ax.vlines(peak[0], peak[1], 1, color="black", linestyle="-", linewidth=0.5)
        ax.text(
            starting_point[0],
            1.1,
            "Starting Point",
            horizontalalignment="center",
            fontsize=fontsize,
        )
        ax.text(
            peak[0],
            1.1,
            "Gradient Peak",
            horizontalalignment="center",
            fontsize=fontsize,
        )

        # 軸の設定
        ax.axis("off")
        ax.set_aspect("equal")
        ax.text(4.5, -0.95, "1 cm", horizontalalignment="center", fontsize=fontsize)
        ax.hlines(-1, 4, 5, color="black", linestyle="-", linewidth=1.5)

        line_collections = [
            LineCollection([], colors="gray", linewidth=1) for _ in range(lines_number)
        ]
        for lc in line_collections:
            ax.add_collection(lc)

        return annotation_boxes, line_collections, starting_point, peak

    annotation_boxes1, line_collections1, starting_point1, peak1 = setup_axis(ax1)
    annotation_boxes2, line_collections2, starting_point2, peak2 = setup_axis(ax2)

    prev_lengths1 = [0] * lines_number
    prev_lengths2 = [0] * lines_number

    prev_min_x1, prev_max_x1, prev_min_y1, prev_max_y1 = None, None, None, None
    prev_min_x2, prev_max_x2, prev_min_y2, prev_max_y2 = None, None, None, None

    total_frames = max(len(lines) for lines in all_lines1) // skip
    time_per_frame = time / total_frames

    # 経過時間を表示するテキスト
    time_text = fig.text(
        0.4,
        0.5,
        "",
        horizontalalignment="center",
        fontsize=fontsize * 1.5,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),
    )

    def animate(frame):
        nonlocal prev_min_x1, prev_max_x1, prev_min_y1, prev_max_y1
        nonlocal prev_min_x2, prev_max_x2, prev_min_y2, prev_max_y2

        def update_lines(
            all_lines,
            line_collections,
            annotation_boxes,
            prev_lengths,
            ax,
            starting_point,
            peak,
            prev_min_x,
            prev_max_x,
            prev_min_y,
            prev_max_y,
        ):
            """ラインとアノテーションボックスの更新"""
            min_x = min(starting_point[0], peak[0])
            max_x = max(starting_point[0], peak[0])
            min_y = min(starting_point[1], peak[1])
            max_y = max(starting_point[1], peak[1])

            for i, lines in enumerate(all_lines):
                prev_length = prev_lengths[i]
                current_length = min(len(lines), frame * skip + 1)
                new_segments = lines[prev_length:current_length]
                if new_segments:
                    segments = list(line_collections[i].get_segments())
                    segments.extend(new_segments)
                    line_collections[i].set_segments(segments)

                if current_length > 0:
                    latest_point = lines[current_length - 1][-1]
                    annotation_boxes[i].xy = (latest_point[0], latest_point[1])

                    if current_length > 1:
                        previous_point = lines[current_length - 2][-1]
                        forward = latest_point[0] >= previous_point[0]
                        annotation_boxes[i].offsetbox = get_imagebox(forward)

                    min_x = min(min_x, latest_point[0])
                    max_x = max(max_x, latest_point[0])
                    min_y = min(min_y, latest_point[1])
                    max_y = max(max_y, latest_point[1])

                prev_lengths[i] = current_length

            # カメラの範囲設定
            target_min_x = min_x - padding
            target_max_x = max_x + padding
            target_min_y = min_y - padding
            target_max_y = max_y + padding

            if prev_min_x is None:
                prev_min_x, prev_max_x, prev_min_y, prev_max_y = (
                    target_min_x,
                    target_max_x,
                    target_min_y,
                    target_max_y,
                )
            else:
                prev_min_x += smooth_factor * (target_min_x - prev_min_x)
                prev_max_x += smooth_factor * (target_max_x - prev_max_x)
                prev_min_y += smooth_factor * (target_min_y - prev_min_y)
                prev_max_y += smooth_factor * (target_max_y - prev_max_y)

            ax.set_xlim(prev_min_x, prev_max_x)
            ax.set_ylim(prev_min_y, prev_max_y)

            return prev_min_x, prev_max_x, prev_min_y, prev_max_y

        prev_min_x1, prev_max_x1, prev_min_y1, prev_max_y1 = update_lines(
            all_lines1,
            line_collections1,
            annotation_boxes1,
            prev_lengths1,
            ax1,
            starting_point1,
            peak1,
            prev_min_x1,
            prev_max_x1,
            prev_min_y1,
            prev_max_y1,
        )
        prev_min_x2, prev_max_x2, prev_min_y2, prev_max_y2 = update_lines(
            all_lines2,
            line_collections2,
            annotation_boxes2,
            prev_lengths2,
            ax2,
            starting_point2,
            peak2,
            prev_min_x2,
            prev_max_x2,
            prev_min_y2,
            prev_max_y2,
        )

        # 必要に応じて濃度マップを更新
        if (
            prev_min_x1 < broad_x_range[0]
            or prev_max_x1 > broad_x_range[1]
            or prev_min_y1 < broad_y_range[0]
            or prev_max_y1 > broad_y_range[1]
        ):
            z = compute_concentration_map(
                [prev_min_x1, prev_max_x1], [prev_min_y1, prev_max_y1], c_mode
            )
            im1.set_data(z)

        if (
            prev_min_x2 < broad_x_range[0]
            or prev_max_x2 > broad_x_range[1]
            or prev_min_y2 < broad_y_range[0]
            or prev_max_y2 > broad_y_range[1]
        ):
            z = compute_concentration_map(
                [prev_min_x2, prev_max_x2], [prev_min_y2, prev_max_y2], c_mode
            )
            im2.set_data(z)

        # フレームに応じた時間の表示
        elapsed_time = frame * time_per_frame
        time_text.set_text(f"Time: {int(elapsed_time)} s")

    # 濃度マップとカラーバーの作成
    cmap = create_light_colormap()
    im1 = ax1.imshow(
        broad_concentration_map,
        cmap=cmap,
        origin="lower",
        extent=[broad_x_range[0], broad_x_range[1], broad_y_range[0], broad_y_range[1]],
        alpha=0.5,
    )

    im2 = ax2.imshow(
        broad_concentration_map,
        cmap=cmap,
        origin="lower",
        extent=[broad_x_range[0], broad_x_range[1], broad_y_range[0], broad_y_range[1]],
        alpha=0.5,
    )

    cbar = fig.colorbar(im1, ax=ax1, shrink=0.6)
    cbar.set_label("Concentration (mM)", fontsize=fontsize)

    cbar = fig.colorbar(im2, ax=ax2, shrink=0.6)
    cbar.set_label("Concentration (mM)", fontsize=fontsize)

    ani = animation.FuncAnimation(
        fig, animate, frames=total_frames, interval=100, repeat=True
    )

    start_time = tm.time()
    ani.save(out_file_path, writer="ffmpeg", dpi=300)
    save_time = tm.time() - start_time
    print(f"Time taken to save the animation: {save_time:.2f} seconds")
