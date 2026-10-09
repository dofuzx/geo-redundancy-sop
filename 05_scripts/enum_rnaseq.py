#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""实验B·第一步：枚举人类CRC RNA-seq系列并探测提交者基因级表达矩阵。

说明：NCBI 统一 counts（*_raw_counts_GRCh38*_NCBI.tsv.gz）仅经 geo/download 端点分发
（bot 防护/reCAPTCHA，无法程序化批量获取），且不在 FTP 目录树内；故与 microarray 臂
保持一致，采用"as-deposited"（提交者自行处理的）基因级矩阵。

流程:
  1) E-utilities esearch (db=gds): colorectal/colon/rectal + Homo sapiens + GSE
  2) esummary 过滤: gdstype == "Expression profiling by high throughput sequencing",
     taxon == "Homo sapiens"（排除人鼠混合系列）
  3) 并发抓取每个系列的 GEO FTP suppl 目录索引, 按 count/tpm/fpkm 关键词挑选基因级矩阵
  4) 输出 rnaseq_enum.tsv (全部 RNA-seq) 与 rnaseq_download_list.tsv (含矩阵的下载清单)

用法:
  python enum_rnaseq.py                 # 全量
  python enum_rnaseq.py --probe-limit 20 --skip-esearch   # 快速探测（调试）
"""
import argparse, concurrent.futures as cf, json, os, re, sys, time
import urllib.parse, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT_DIR = r"D:/AIwork/GEO_RNAseq"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
RNA_TYPE = "Expression profiling by high throughput sequencing"
UA = {"User-Agent": "geo-redundancy-sop/1.0 (experiment-B rnaseq enumeration)"}


def http_json(url, data=None, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=UA)
            return json.loads(urllib.request.urlopen(req, timeout=60).read())
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def esearch(term, retmax=12000):
    q = urllib.parse.quote(term)
    url = f"{EUTILS}/esearch.fcgi?db=gds&term={q}&retmax={retmax}&retmode=json"
    d = http_json(url)
    n = int(d["esearchresult"]["count"])
    ids = d["esearchresult"]["idlist"]
    if n > len(ids):
        print(f"  [warn] count={n} > retmax={retmax}, 截断")
    return ids


def esummary_all(ids, batch=500):
    out = {}
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        data = urllib.parse.urlencode({"db": "gds", "retmode": "json",
                                       "id": ",".join(chunk)}).encode()
        d = http_json(f"{EUTILS}/esummary.fcgi", data=data)
        res = d["result"]
        for uid in res["uids"]:
            out[uid] = res[uid]
        print(f"  esummary {min(i+batch, len(ids))}/{len(ids)}", flush=True)
    return out


def gse_bucket(gse_num):
    """GEO FTP bucket: GSE20916 -> GSE20nnn"""
    s = str(gse_num)
    return "GSE" + (s[:-3] if len(s) >= 4 else "0") + "nnn"


# 接受的提交者基因级矩阵: 名字含 count/tpm/fpkm, 非转录本级, 常规表格扩展名
CAND_RE = re.compile(r'(count|tpm|fpkm)', re.I)
BAD_RE = re.compile(r'(transcript|mirna|mir_|isoform|_sj_|junc|splice|peak|methylation)', re.I)
EXTS = ('.txt.gz', '.csv.gz', '.tsv.gz', '.txt', '.csv', '.tsv', '.xz', '.gz')
PREF = [('raw', 'count'), ('count',), ('tpm',), ('fpkm',)]


def pick_candidate(names):
    cands = [n for n in names
             if CAND_RE.search(n) and not BAD_RE.search(n)
             and n.lower().endswith(EXTS) and not n.lower().endswith('.tar')]
    if not cands:
        return None, []
    def rank(n):
        low = n.lower()
        for i, keys in enumerate(PREF):
            if all(k in low for k in keys):
                return i
        return len(PREF)
    cands.sort(key=lambda n: (rank(n), -len(n)))
    return cands[0], cands


def probe_suppl(gse):
    """返回 (gse, chosen_file or None, candidates_joined, error or '')"""
    url = f"{FTP}/{gse_bucket(int(gse[3:]))}/{gse}/suppl/"
    last = ""
    for i in range(3):
        try:
            req = urllib.request.Request(url, headers=UA)
            html = urllib.request.urlopen(req, timeout=45).read().decode("utf-8", "replace")
            names = [n for n in re.findall(r'href="([^"]+)"', html)
                     if not n.startswith('/') and not n.startswith('http')]
            chosen, cands = pick_candidate(names)
            return (gse, chosen, ";".join(cands), "")
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
            time.sleep(1.5 * (i + 1))
    return (gse, None, "", last)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-limit", type=int, default=0, help="只探测前N个（调试用）")
    ap.add_argument("--skip-esearch", action="store_true",
                    help="跳过esearch/esummary，直接用已有 rnaseq_enum.tsv 探测suppl")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    if not args.skip_esearch:
        print("[1/3] esearch ...", flush=True)
        term = ('(colorectal[All Fields] OR "colon cancer"[All Fields] '
                'OR "rectal cancer"[All Fields]) AND Homo sapiens[Organism] AND gse[ETYP]')
        ids = esearch(term)
        print(f"  esearch 命中 {len(ids)} 条", flush=True)
        print("[2/3] esummary + 过滤 RNA-seq / human ...", flush=True)
        docs = esummary_all(ids)
        rows = []
        for uid, d in docs.items():
            if d.get("entrytype") != "GSE":
                continue
            if d.get("gdstype") != RNA_TYPE:
                continue
            tax = (d.get("taxon") or "").strip()
            if tax != "Homo sapiens":
                continue
            n = d.get("n_samples") or 0
            try:
                n = int(n)
            except (TypeError, ValueError):
                n = 0
            rows.append({"gse": d.get("accession"), "n_samples": n,
                         "title": (d.get("title") or "")[:200],
                         "pdat": d.get("pdat") or "",
                         "suppfile": d.get("suppfile") or ""})
        rows.sort(key=lambda r: r["gse"])
        with open(os.path.join(OUT_DIR, "rnaseq_enum.tsv"), "w", encoding="utf-8", newline="") as f:
            w = __import__("csv").DictWriter(f, fieldnames=["gse", "n_samples", "title", "pdat", "suppfile"], delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"  人类CRC RNA-seq 系列: {len(rows)} 个 -> rnaseq_enum.tsv", flush=True)
    else:
        with open(os.path.join(OUT_DIR, "rnaseq_enum.tsv"), encoding="utf-8") as f:
            rows = list(__import__("csv").DictReader(f, delimiter="\t"))
        print(f"  载入已有清单 {len(rows)} 个", flush=True)

    todo = [r["gse"] for r in rows]
    if args.probe_limit:
        todo = todo[:args.probe_limit]
    print(f"[3/3] 探测 {len(todo)} 个系列的 suppl 目录 (NCBI统一counts) ...", flush=True)
    found, miss, err = [], 0, 0
    cand_map = {}
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=16) as ex:
        futs = {ex.submit(probe_suppl, g): g for g in todo}
        for k, fut in enumerate(cf.as_completed(futs), 1):
            gse, fname, cands, e = fut.result()
            if e:
                err += 1
            if fname:
                found.append((gse, fname))
                cand_map[gse] = cands
            else:
                miss += 1
            if k % 100 == 0 or k == len(todo):
                print(f"  {k}/{len(todo)}  有矩阵={len(found)}  无={miss}  错误={err}  "
                      f"({time.time()-t0:.0f}s)", flush=True)
    with open(os.path.join(OUT_DIR, "rnaseq_with_counts.tsv"), "w", encoding="utf-8", newline="") as f:
        w = __import__("csv").writer(f, delimiter="\t")
        w.writerow(["gse", "counts_file", "candidates"])
        for gse, fname in sorted(found):
            w.writerow([gse, fname, cand_map.get(gse, "")])
    fmap = dict(found)
    n = [r for r in rows if r["gse"] in fmap]
    with open(os.path.join(OUT_DIR, "rnaseq_download_list.tsv"), "w", encoding="utf-8", newline="") as f:
        w = __import__("csv").DictWriter(f, fieldnames=["gse", "n_samples", "title", "counts_file"],
                                         delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(n, key=lambda r: -int(r["n_samples"] or 0)):
            w.writerow({**r, "counts_file": fmap[r["gse"]]})
    print(f"[完成] 含提交者基因级矩阵的系列: {len(found)} 个 -> rnaseq_download_list.tsv", flush=True)


if __name__ == "__main__":
    main()
