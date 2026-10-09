# Tier4（领域级 200+）数据下载操作指南

> 目标：把你已有的 61 系列（Tier1+2+3）扩展到**领域级 250–300 系列**，用于估算 CRCRedundancy 的真实广度。  
> 适用：本机可访问 GEO（`ftp.ncbi.nlm.nih.gov`）的环境。生成时间 2026-10-09。

---

## 零、先看清全局：已生成的两个清单

| 清单文件                         | 条目      | 平台           | 说明                | 建议         |
| ---------------------------- | ------- | ------------ | ----------------- | ---------- |
| `crc_gse_tier4_highconf.txt` | **109** | GPL570/GPL96 | 标题明确是 CRC/结肠/直肠队列 | ✅ **先下这批** |
| `crc_gse_tier4_full.txt`     | **247** | GPL570/GPL96 | 含泛组织/细胞系/血清大队列    | 第二批，或按需    |

**为什么分两批**：`--full` 里最大的几个其实不是 CRC 专属队列，而是碰巧含 CRC 样本的泛队列：

| GSE       | 样本   | 标题                                                   | 风险                 |
| --------- | ---- | ---------------------------------------------------- | ------------------ |
| GSE203024 | 2845 | Gene Expression Data from Human **Peripheral Blood** | 外周血，非肿瘤组织，会稀释冗余率分母 |
| GSE2109   | 2158 | Expression Project for Oncology (expO)               | 泛癌多组织混合            |
| GSE7307   | 677  | Human body index                                     | 多器官普查              |
| GSE57083  | 627  | AstraZeneca **cell lines**                           | 细胞系，非患者样本          |

高置信清单里最大的是 `GSE14333`（290 例原发 CRC）、`GSE161158`(250)、`GSE101896`(190)、`GSE14095`(189)、`GSE13294`(155) ——这些才真正扩充领域语料。

### 📊 实测下载规划（已逐条 HEAD 探测，2026-10-09）

| 清单                           | 条目  | 可达性               | 总下载量(压缩)       | 平均/个     |
| ---------------------------- | --- | ----------------- | -------------- | -------- |
| `crc_gse_tier4_highconf.txt` | 109 | 104/109 探测到（其余见下） | **≈ 0.97 GiB** | 9.5 MiB  |
| `crc_gse_tier4_full.txt`     | 247 | **247/247 全部可达**  | **≈ 3.46 GiB** | 14.7 MiB |

- 体积远小于预期，**高置信首批不到 1 GiB**，很快就下完。
- 单个最大的是 `GSE203024`（**515.6 MiB**，外周血泛队列）；高置信里最大的是 `GSE14333`（80.8 MiB）、`GSE161158`(69.9 MiB)。
- 注：初次并发探测时 `GSE30292/GSE13059/GSE24795/GSE54483/GSE7678/GSE22061/GSE31084/GSE16080` 报失败，复查确认是**并发瞬时网络错误**，重新请求全部可达（合计仅 34 MiB）。下载器本身有重试机制，不影响。

---

## 一、三步走（复制即用）

### 第 1 步：先干跑，确认 URL 与路径（零下载）

"D:\AIwork\GEO\Tire4"

期望输出：`已解析 109 个` + 逐条 `桶=GSE__nnn OK`。若某条 FAIL，说明该 GSE 的 matrix 路径异常，跳过即可。

### 第 2 步：正式下载（推荐 3 并发，礼貌 NCBI）

```bash
python download_tier3_4.py --list crc_gse_tier4_highconf.txt --out-root D:/AIwork/GEO_T4 --workers 3
```

运行界面（已带实时进度）：

"D:\AIwork\GEO\Tire4"

### 第 3 步：下完第二批改同一个根目录（会自动跳过已存在的）

```bash
python download_tier3_4.py --list crc_gse_tier4_full.txt --out-root D:/AIwork/GEO_T4 --workers 3
```

产物：`download_manifest.json`（成功清单）、`download_failed.txt`（失败清单，可直接当 `--list` 重下）。

---

## 二、下载后的目录结构

```
D:/AIwork/GEO_T4/
├── GSE14333/GSE14333_series_matrix.txt.gz
├── GSE161158/GSE161158_series_matrix.txt.gz
├── ...
├── download_manifest.json
└── download_failed.txt
```

**重要**：请保证新目录与已有 Tier3 目录（`D:/AIwork/GEO`）分开，方便分层复算。

---

## 三、下载后交给我处理的衔接命令（供参考）

```bash
PY=C:/Users/fuzex/.workbuddy/binaries/python/envs/default/Scripts/python.exe
cd C:/Users/fuzex/WorkBuddy/2026-10-07-16-35-27

# 1) 完整性校验（必须先做！上次 Tier3 有 2 个文件下载截断，导致转换中途崩溃）
"$PY" check_gz_integrity.py            # 需把根目录改为 D:/AIwork/GEO_T4

# 2) 流式转换成 SOP 格式（零依赖，不会 OOM）
"$PY" geo_matrix_to_sop_stream.py --root D:/AIwork/GEO_T4 --out D:/AIwork/GEO_T4_sop

# 3) 剔除 0 基因的 RNA-seq 系列 + SuperSeries 后，与 Tier1+2+3 合并检测
"$PY" analyze_tier.py \
    --indirs D:/AIwork/GEO_T4_sop D:/AIwork/GEO_sop_v2 E:/CRC_GEO_2026T1_sop_v2 E:/CRC_GEO_2026T2_sop_v2 \
    --out D:/AIwork/GEO_ALL_dedup --name "Field-wide" --pos1 GSE41258 --pos2 GSE68468
```

---

## 四、重新/增量生成清单（清单会过期，可随时刷新）

```bash
# 领域级 250 条（默认平台 GPL570 + GPL96）
python gen_tier4_list.py --target 250

# 只要高置信（标题必须命中 CRC 关键词）
python gen_tier4_list.py --target 300 --require-title

# 换个平台组合，或降低样本门槛
python gen_tier4_list.py --platforms 570 96 97 --min-samples 10 --target 300
```

`gen_tier4_list.py` 通过 NCBI E-utilities **程序化**枚举，**不要手工抄 200+ 个编号**（极易出错且不可复现）。  
已内置：按系列自动剔除既用编号、自动排除 SuperSeries、NCBI 限速（无 API key 时 3 req/s）与失败重试。

有 NCBI API key 可提速 10 倍：

```bash
python gen_tier4_list.py --target 250 --api-key 你的KEY
```

---

## 五、注意事项

1. **并发别开太高**：NCBI 对高频请求会封 IP，`--workers` 保持 ≤3。
2. **务必先干跑**：`--dry-run` 能在零下载的情况下验证 247 条路径是否正确。
3. **下载完必须校验 gzip 完整性**：上次 Tier3 有 `GSE22598`/`GSE32323` 被截断在 7MB 却仍被记为成功，直接转换会中途崩溃。
4. **再来一遍拦截器**：体积大的系列（如 xenograft panel >500MB）耗时可能数十分钟，进度条会显示内部进度，耐心等待，不要中断（中断会产生截断文件）。
5. 若某条反复失败，把编号加进 `--exclude`，或直接删掉该行。
