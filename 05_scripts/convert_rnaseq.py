#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""实验B·第三步：异构提交者矩阵 → 统一 SOP 格式 (genes × samples CSV + meta TSV)。

鲁棒解析: gzip/plain, 分隔符嗅探(\\t/,/;), 引号清理, 样本名去重;
校验: 样本列 4-300, 基因行 >=5000, 基因ID唯一率>=0.9;
产出: D:/AIwork/GEO_RNAseq/sop/{GSE}.csv 与 {GSE}_meta.tsv (platform=RNA-seq_count)
"""
import csv, glob, gzip, io, os, re, sys

SRC = r"D:/AIwork/GEO_RNAseq/counts"
DST = r"D:/AIwork/GEO_RNAseq/sop"
os.makedirs(DST, exist_ok=True)


def sniff_delim(line):
    for d in ("\t", ",", ";"):
        if line.count(d) >= 2:
            return d
    return None          # 空白分隔


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if path.endswith(".xz"):
        import lzma
        return lzma.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def fnum(v):
    try:
        if v in ("", "NA", "NaN", "nan", "null", "NULL", "-"):
            return ""
        return round(float(v), 4)
    except ValueError:
        return ""


# 注释列黑名单（归一化后精确匹配）
_ANORM = lambda s: re.sub(r"[^a-z0-9]", "", s.lower())
ANNOT_NAMES = {
    "chrom", "chromosome", "chr", "chromname", "chrname", "start", "end", "strand",
    "length", "lengthgene", "gc", "gccontent", "genename", "genesymbol", "symbol",
    "biotype", "genebiotype", "genetype", "type", "description", "genedescription",
    "tffamily", "externalgenename", "entrez", "entrezid", "entrezgene", "entrezgeneid",
    "refseq", "refseqid", "refseqmrna", "band", "cytoband", "genestart", "geneend",
    "genestrand", "genelength", "genchr", "genestartbp", "geneendbp", "probeid",
    "meanexpr", "medianexpr", "log2fc", "logfc", "padj", "pvalue", "fdr", "pval",
    "qvalue", "featuretype", "seqname", "position",
    # 补充: featureCounts/HtSeq 等常见注释列变体
    "exonicgenesize", "exonicgenesizes", "transcriptlength", "transcriptlengths",
    "pos", "stop", "width", "geneid", "id", "name", "locus", "contig", "coords",
    "genomiclocation", "exonlength", "efflength", "bases", "basesmapped",
}
LOG_MARKERS = ("command:", "program:", "version:", "options:", "started", "runinfo",
               "member", "alignment", "output:", "input:")


def looks_like_header(line, d):
    """判断是否为真表头: 含 >=5 个非空字段且不含日志标记。"""
    parts = [p.strip() for p in (line.split(d) if d else line.split())]
    parts = [p for p in parts if p]
    if len(parts) < 5:
        return False
    low = " ".join(parts).lower()
    return not any(m in low for m in LOG_MARKERS)


def find_header(f, d):
    """跳过前导 '#' 注释行与日志行, 返回 (header_line, consumed)。"""
    consumed = 0
    while True:
        pos = f.tell()
        line = f.readline()
        if not line:
            f.seek(pos)
            return None, consumed
        s = line.strip()
        consumed += 1
        if not s or s.startswith("#"):
            continue
        if looks_like_header(s, d):
            f.seek(pos)          # main loop 从该行重新读
            return s, consumed
        # 日志行 -> 继续向下找
    return None, consumed


def pick_sample_cols(path, d, hdr):
    """两遍法确定真正的样本列: 表头黑名单 + 前3000数据行的低基数检测。
    返回 (keep_idx 相对第1列起的列索引, n_dropped)。"""
    cand = list(range(1, len(hdr)))          # 绝对列索引
    keep = [i for i in cand if _ANORM(hdr[i]) not in ANNOT_NAMES]
    # 低基数检测: 注释列(染色体/链/生物型等)取值集合极小
    if keep:
        lim = 3000
        seen = {i: set() for i in keep}
        with open_maybe_gz(path) as f:
            first = f.readline()
            dd = sniff_delim(first)
            f.seek(0)
            find_header(f, dd)
            for k, line in enumerate(f):
                if k >= lim:
                    break
                if not line.strip():
                    continue
                parts = line.rstrip("\r\n").split(d) if d else line.split()
                for i in keep:
                    if i < len(parts):
                        v = parts[i].strip()
                        if v != "":
                            seen[i].add(v)
        keep = [i for i in keep if len(seen[i]) > 200]
    return keep, len(cand) - len(keep)


def convert(path, gse):
    # 二进制格式守卫: .xls(OLE2)/.xlsx(ZIP) 冒充文本扩展名（.gz 需看解压后字节）
    try:
        if path.endswith(".gz"):
            with gzip.open(path, "rb") as fb:
                magic = fb.read(8)
        else:
            with open(path, "rb") as fb:
                magic = fb.read(8)
        if magic[:4] == b"\xd0\xcf\x11\xe0" or magic[:4] == b"PK\x03\x04":
            return "binary_excel_not_text"
    except OSError:
        return "corrupt_gzip"
    with open_maybe_gz(path) as f:
        first = f.readline().rstrip("\r\n")
        if not first.strip():
            return "empty_header"
        d = sniff_delim(first)
        f.seek(0)
        raw_header, _skipped = find_header(f, d)
        if raw_header is None:
            return "no_valid_header"
        header = raw_header.rstrip("\r\n")
        d = sniff_delim(header)
        hdr = [c.strip().strip('"') for c in (header.split(d) if d else header.split())]
        ncol = len(hdr)
        if ncol < 5:                      # 首列基因 + >=4 样本(含可能注释列)
            return f"too_few_columns({ncol})"
        keep_idx, n_annot = pick_sample_cols(path, d, hdr)
        samples = [hdr[i] for i in keep_idx]
        if len(samples) < 4:
            return f"too_few_real_samples({len(samples)}, annot_cols={n_annot})"
        # 样本名清洗与去重
        seen = {}
        clean = []
        for s in samples:
            s = re.sub(r"[\s\"']+", "_", s)[:80] or "c"
            if s in seen:
                seen[s] += 1
                s = f"{s}.{seen[s]}"
            else:
                seen[s] = 0
            clean.append(s)
        if len(clean) > 300:
            return f"sample_count_out_of_range({len(clean)})"
        out_csv = os.path.join(DST, gse + ".csv")
        n_gene = 0
        ids = set()
        dup = 0
        with open(out_csv, "w", newline="", encoding="utf-8") as fo:
            w = csv.writer(fo)
            w.writerow(["gene_id"] + clean)
            for ln, line in enumerate(f):
                if not line.strip():
                    continue
                parts = line.rstrip("\r\n").split(d) if d else line.split()
                if len(parts) < max(keep_idx) + 1:
                    parts += [""] * (max(keep_idx) + 1 - len(parts))
                gid = parts[0].strip().strip('"')
                if not gid:
                    continue
                if gid in ids:
                    dup += 1
                    continue                     # 重复基因行: 保留首个
                ids.add(gid)
                # 少量行带注释行(以#开头)——跳过
                if ln > 0 and gid.startswith("#"):
                    ids.discard(gid)
                    continue
                vals = [fnum(parts[i].strip()) for i in keep_idx]
                if all(v == "" for v in vals):
                    continue
                w.writerow([gid] + vals)
                n_gene += 1
                if n_gene > 200000:
                    fo.close(); os.remove(out_csv)
                    return "too_many_rows(>200k, 疑似单细胞)"
        uniq = len(ids) / max(1, n_gene)
        if n_gene < 5000:
            os.remove(out_csv)
            return f"too_few_genes({n_gene})"
        if uniq < 0.9:
            os.remove(out_csv)
            return f"gene_ids_not_unique({uniq:.2f})"
    # meta.tsv
    with open(os.path.join(DST, gse + "_meta.tsv"), "w", newline="", encoding="utf-8") as fm:
        w = csv.writer(fm, delimiter="\t")
        w.writerow(["sample_id", "gsm", "platform"])
        for s in clean:
            w.writerow([s, s, "RNA-seq_count"])
    return f"ok({n_gene} genes x {len(clean)} samples, dup={dup}, annot_cols_dropped={n_annot})"


def main():
    files = sorted(glob.glob(os.path.join(SRC, "*__*")))
    print(f"待转换 {len(files)} 个矩阵", flush=True)
    okn = 0
    with open(os.path.join(DST, "convert_report.tsv"), "w", encoding="utf-8", newline="") as fr:
        w = csv.writer(fr, delimiter="\t")
        w.writerow(["gse", "file", "result"])
        for k, path in enumerate(files, 1):
            base = os.path.basename(path)
            gse = base.split("__", 1)[0]
            try:
                res = convert(path, gse)
            except Exception as e:
                res = f"EXC:{type(e).__name__}:{str(e)[:80]}"
            if res.startswith("ok"):
                okn += 1
            w.writerow([gse, base, res])
            if k % 25 == 0 or k == len(files):
                print(f"  {k}/{len(files)}  ok={okn}  last={gse}:{res[:60]}", flush=True)
    print(f"[完成] 转换成功 {okn}/{len(files)} -> {DST}", flush=True)


if __name__ == "__main__":
    main()
