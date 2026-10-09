#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""扫描 D:/AIwork/GEO 下所有 <GSE>_series_matrix.txt.gz 的 gzip 完整性。
输出损坏(截断)文件列表，便于用户定向重下。"""
import gzip, os, sys

ROOT = "D:/AIwork/GEO"

def main():
    bad = []
    good = 0
    for name in sorted(os.listdir(ROOT)):
        d = os.path.join(ROOT, name)
        if not (os.path.isdir(d) and name.startswith("GSE")):
            continue
        gz = os.path.join(d, name + "_series_matrix.txt.gz")
        if not os.path.exists(gz):
            continue
        size = os.path.getsize(gz)
        try:
            with gzip.open(gz, "rt", encoding="utf-8", errors="replace") as f:
                n = 0
                for _ in f:
                    n += 1
                good += 1
                print(f"[OK]   {name:14s} {size/1024/1024:7.1f}MB  {n} 行")
        except Exception as e:
            bad.append((name, size, str(e)))
            print(f"[坏]   {name:14s} {size/1024/1024:7.1f}MB  -> {e}")
    print("\n==== 汇总 ====")
    print(f"完好: {good}  损坏: {len(bad)}")
    if bad:
        print("损坏清单(需重下):")
        for name, size, e in bad:
            print(f"  {name}  ({size/1024/1024:.1f}MB)  {e}")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
