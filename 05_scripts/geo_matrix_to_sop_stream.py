#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEO series_matrix.txt.gz -> SOP 格式【流式】转换器

与技能版 geo_matrix_to_sop.py 输出完全兼容：
  <NAME>.csv       表达矩阵：首列 ID_REF，其余列样本(GSM)，行=基因/探针
  <NAME>_meta.tsv  样本元数据：tab 分隔，含 sample_id / gsm / series_title / platform / 各 !Sample_* 字段

【为何需要流式】技能版把整个表达表 accumulate 到 list 里再写盘，
遇到巨型系列（如 GSE39582，165MB 压缩 / 54675 探针 × 数百样本）会 OOM/被杀。
本版边解析边写盘，内存只占：元数据 + 当前一行。

【命名规则】
  单个 matrix 映射到某 GSE        -> 用原 GSE 名（如 GSE17536）
  同一 GSE 有多个 matrix(多平台)  -> 派生平台后缀名（如 GSE35982_GPL14767、GSE17538_GPL570）
  这样既复刻既有 Tier2 双平台命名，又能正确处理 Tier1 的 GSE17538(GPL570+GPL1261)。

依赖：仅标准库。Windows 路径请用 E:/ 形式。
"""
import argparse
import csv
import gzip
import os
import re
import sys


def _unquote(s):
    s = s.strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    return s


def derived_name(gz_path):
    """GSE17538-GPL570_series_matrix.txt.gz -> GSE17538_GPL570
       GSE110224_series_matrix.txt.gz        -> GSE110224"""
    base = os.path.basename(gz_path)
    name = re.sub(r"[._\-]?series_matrix\.txt\.gz$", "", base, flags=re.I)
    name = name.replace("-", "_").strip()
    return name


def collect_targets(root):
    """递归找出所有 series_matrix，按 GSE 归组产出 (path, out_name) 列表"""
    files = []
    for cur, _dirs, fs in os.walk(root):
        for fn in fs:
            if fn.lower().endswith("series_matrix.txt.gz"):
                files.append(os.path.join(cur, fn))
    files = sorted(set(files))

    groups = {}
    for p in files:
        m = re.search(r"(GSE\d+)", os.path.basename(p))
        if not m:
            continue
        groups.setdefault(m.group(1), []).append(p)

    targets = []
    for gse, paths in groups.items():
        # 同 GSE 内按体积降序：主数据集(通常更大)优先，输出顺序稳定
        paths.sort(key=lambda p: os.path.getsize(p), reverse=True)
        if len(paths) == 1:
            targets.append((paths[0], gse))
        else:
            for p in paths:
                targets.append((p, derived_name(p)))
    return targets


def process_one(gz_path, out_dir, name):
    """流式转换单个 series_matrix -> <name>.csv + <name>_meta.tsv"""
    meta = {}           # field -> [values per sample]
    meta_order = []     # 字段顺序
    series_title = ""
    platform = ""

    csv_path = os.path.join(out_dir, name + ".csv")
    tsv_path = os.path.join(out_dir, name + "_meta.tsv")

    header = None
    samples = []
    nrows = 0
    saw_table = False

    os.makedirs(out_dir, exist_ok=True)
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f, \
            open(csv_path, "w", newline="", encoding="utf-8") as out:
        w = csv.writer(out)
        in_table = False

        for raw in f:
            line = raw.rstrip("\n")

            if in_table:
                if header is None:
                    header = line.split("\t")
                    samples = header[1:]
                    w.writerow(header)
                    continue
                if line.startswith("!series_matrix_table_end"):
                    in_table = False
                    continue
                parts = line.split("\t")
                if not parts or parts == [""]:
                    continue
                # 对齐列数，避免 ragged CSV 破坏下游 numpy 加载
                if len(parts) > len(header):
                    parts = parts[:len(header)]
                elif len(parts) < len(header):
                    parts = parts + [""] * (len(header) - len(parts))
                w.writerow(parts)
                nrows += 1
                if nrows % 20000 == 0:
                    print(f"    · {name} 已写 {nrows} 行", flush=True)
                continue

            # ---- 表之前的元数据区 ----
            if line.startswith("!Series_title"):
                series_title = _unquote(line.split("\t", 1)[1]) if "\t" in line else ""
            elif line.startswith("!Series_platform_id"):
                platform = _unquote(line.split("\t", 1)[1]) if "\t" in line else ""
            elif line.startswith("!Sample_"):
                field = line[1:].split("\t", 1)[0]
                vals = line.split("\t")[1:] if "\t" in line else []
                meta[field] = [_unquote(v) for v in vals]
                if field not in meta_order:
                    meta_order.append(field)
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                saw_table = True
                header = None

    if not saw_table:
        # 完全没有表达表：移除空壳，标记跳过（RNA-seq 只给 SRA 链接时常见）
        try:
            os.remove(csv_path)
        except OSError:
            pass
        print(f"[跳过] {name}: 未找到表达表", flush=True)
        return None

    # 优先使用逐样本的 !Sample_platform_id，因其与实际探针集一致：
    #   系列级 !Series_platform_id 会误标（如 GSE17538 系列级=GPL1261，
    #   但其 GPL570 子集的实际平台为 GPL570；GSE35982 双平台同理）。
    splat = meta.get("Sample_platform_id", [])
    cols = ["sample_id", "gsm", "series_title", "platform"] + meta_order
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for i, sid in enumerate(samples):
            p_i = splat[i] if (i < len(splat) and splat[i]) else platform
            row = [sid, sid, series_title, p_i]
            for fld in meta_order:
                vlist = meta.get(fld, [])
                row.append(vlist[i] if i < len(vlist) else "")
            w.writerow(row)

    print(f"[OK] {name}: {len(samples)} 样本 × {nrows} 基因  platform={platform}", flush=True)
    return name


def main(argv=None):
    ap = argparse.ArgumentParser(description="GEO series_matrix.gz -> SOP 格式（流式，支持超大矩阵）")
    ap.add_argument("--root", help="根目录，含 GSExxxx 子目录（可递归含 processed/ 等层级）")
    ap.add_argument("--file", help="直接指定单个 *_series_matrix.txt.gz")
    ap.add_argument("--out", help="输出扁平目录")
    args = ap.parse_args(argv)

    if args.file:
        targets = [(args.file, derived_name(args.file))]
        out_dir = args.out or "./sop_out"
    elif args.root:
        targets = collect_targets(args.root)
        out_dir = args.out or (os.path.dirname(os.path.abspath(args.root)) + "/" +
                               os.path.basename(args.root.rstrip("/\\")) + "_sop")
    else:
        sys.exit("[错误] 需提供 --root 或 --file")

    print(f"发现 {len(targets)} 个系列；输出目录: {out_dir}", flush=True)
    os.makedirs(out_dir, exist_ok=True)

    done = 0
    for i, (p, name) in enumerate(targets, 1):
        sz = os.path.getsize(p) / (1024 * 1024)
        print(f"[{i}/{len(targets)}] {name}  ({sz:.1f} MB gz)", flush=True)
        try:
            if process_one(p, out_dir, name):
                done += 1
        except Exception as e:
            print(f"[FAIL] {name}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
    print(f"[完成] 成功转换 {done}/{len(targets)} -> {out_dir}", flush=True)


if __name__ == "__main__":
    main()
