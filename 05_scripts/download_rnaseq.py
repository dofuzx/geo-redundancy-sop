#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""实验B·第二步：下载 RNA-seq 基因级表达矩阵（提交者提供，FTP suppl）。

读取 enum_rnaseq.py 产出的 rnaseq_download_list.tsv，并发下载到
D:/AIwork/GEO_RNAseq/counts/{GSE}__{文件名}，支持断点续传与大小上限。
"""
import argparse, concurrent.futures as cf, csv, os, sys, time
import urllib.request

OUT_DIR = r"D:/AIwork/GEO_RNAseq/counts"
LIST_TSV = r"D:/AIwork/GEO_RNAseq/rnaseq_download_list.tsv"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
MAX_BYTES = 300 * 1024 * 1024      # 单文件上限 300MB（巨型单细胞矩阵直接跳过）
UA = {"User-Agent": "geo-redundancy-sop/1.0 (experiment-B rnaseq download)"}


def gse_bucket(gse_num):
    s = str(gse_num)
    return "GSE" + (s[:-3] if len(s) >= 4 else "0") + "nnn"


def download_one(gse, fname):
    dest = os.path.join(OUT_DIR, f"{gse}__{fname}")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return (gse, "cached", os.path.getsize(dest))
    part = dest + ".part"
    url = f"{FTP}/{gse_bucket(gse[3:])}/{gse}/suppl/{urllib.request.quote(fname)}"
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=180) as r, open(part, "wb") as f:
                total = 0
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
                    total += len(b)
                    if total > MAX_BYTES:
                        f.close(); os.remove(part)
                        return (gse, "too_big", total)
            if total < 1000:
                os.remove(part)
                return (gse, "tiny_file", total)
            os.replace(part, dest)
            return (gse, "ok", total)
        except Exception as e:
            if os.path.exists(part):
                os.remove(part)
            if attempt == 3:
                return (gse, f"fail:{type(e).__name__}", 0)
            time.sleep(2 * (attempt + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="只下载前N个（按样本数降序）")
    ap.add_argument("--min-samples", type=int, default=4)
    ap.add_argument("--max-samples", type=int, default=300, help="排除单细胞级超大队列")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(LIST_TSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rows = [r for r in rows if args.min_samples <= int(r["n_samples"] or 0) <= args.max_samples]
    if args.limit:
        rows = rows[:args.limit]
    print(f"待下载 {len(rows)} 个系列 (样本数 {args.min_samples}-{args.max_samples})", flush=True)
    ok = cached = fail = 0
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(download_one, r["gse"], r["counts_file"]): r for r in rows}
        for k, fut in enumerate(cf.as_completed(futs), 1):
            gse, st, size = fut.result()
            if st in ("ok", "cached"):
                ok += 1
                if st == "cached":
                    cached += 1
            else:
                fail += 1
                print(f"  [{st}] {gse} size={size}", flush=True)
            if k % 20 == 0 or k == len(rows):
                print(f"  {k}/{len(rows)}  ok={ok}(缓存{cached})  fail={fail}  "
                      f"{(time.time()-t0)/60:.1f}min", flush=True)
    print(f"[完成] ok={ok} fail={fail} -> {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
