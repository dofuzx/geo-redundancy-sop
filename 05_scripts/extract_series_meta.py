#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从源 series_matrix.txt.gz 抽取【系列级】元数据（所有 !Series_* 字段），
供 SOP step-4 语义复核使用。

meta.tsv 里只有 !Sample_* + Series_title/platform，缺了做语义复核最关键的信息
（summary / overall_design / submitter / pubmed_id / 机构），这些都在 gz 的表头里。
本脚本只读到 !series_matrix_table_begin 就停，因此即使面对 GSE39582(165MB) 也秒级完成。
"""
import gzip
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from geo_matrix_to_sop_stream import collect_targets  # noqa: E402

ROOTS = {
    "Tier3": "D:/AIwork/GEO",
    "Tier1": "E:/CRC_GEO_2026T1",
    "Tier2": "E:/CRC_GEO_2026T2",
}


def _clean(v):
    v = v.strip()
    if len(v) >= 2 and v[0] == '"' and v[-1] == '"':
        v = v[1:-1]
    return v


def series_meta(gz):
    acc = defaultdict(list)
    with gzip.open(gz, "rt", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Series_"):
                continue
            parts = line[1:].split("\t")
            k = parts[0]
            vals = [_clean(v) for v in parts[1:]]
            vals = [v for v in vals if v]
            if vals:
                acc[k].extend(vals)
    return {k: (v[0] if len(v) == 1 else v) for k, v in acc.items()}


def main():
    index = {}
    for tier, root in ROOTS.items():
        if not os.path.isdir(root):
            print(f"[警告] 根目录不存在: {root}", file=sys.stderr)
            continue
        for p, name in collect_targets(root):
            index.setdefault(name, {"tier": tier, "gz": p})

    meta = {}
    for name, info in sorted(index.items()):
        try:
            m = series_meta(info["gz"])
            m["_source"] = {"tier": info["tier"], "gz": info["gz"]}
            meta[name] = m
        except Exception as e:
            meta[name] = {"_error": f"{type(e).__name__}: {e}",
                          "_source": {"tier": info["tier"], "gz": info["gz"]}}

    out = os.path.join(HERE, "series_metadata.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=1)
    print(f"[完成] 已抽取 {len(meta)} 个系列的 series 级元数据 -> {out}")

    # 打印关键字段速览
    keys = ["Series_title", "Series_pubmed_id", "Series_submission_date",
            "Series_last_update_date"]
    print("\n=== 关键字段速览 ===")
    for name, m in sorted(meta.items()):
        if "_error" in m:
            print(f"  {name}: ERROR {m['_error']}")
            continue
        title = str(m.get("Series_title", ""))[:70]
        pm = m.get("Series_pubmed_id", "")
        pm = (",".join(pm) if isinstance(pm, list) else pm)
        sub = m.get("Series_submission_date", "")
        contrib = m.get("Series_contributor_1", "")
        print(f"  {name:22s} PMID={pm:12s} 提交={str(sub):10s} | {title}")
        print(f"  {'':22s} 主要贡献者={contrib}")


if __name__ == "__main__":
    main()
