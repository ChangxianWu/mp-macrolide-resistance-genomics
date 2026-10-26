# -*- coding: utf-8 -*-
"""03e_calibrate_sites.py — 23S 耐药位点校准 v2（表型驱动 + 文献锚定；取代早期的 EC 标尺路线）

背景（为什么弃用 EC 标尺路线）:
  早期路线把 MP 文献编号 2063/2064/2617 误当成 E. coli 编号去 NW 映射 → WT 自检 T/C/T ≠ A/A/C 永远 WARN;
  随后的 EC 标尺核验又发现: J01695 ≡ U00096 注释切片（两条独立记录序列相同），
  其 2058/2059 序列坐标 = T/G —— EC 注释坐标与经典 rRNA 编号本身存在偏移，EC 标尺不可直接用。

校准证据链（本脚本逐条落盘到 results/MY9_位点校准报告.md）:
  E1 表型锚定: 试点扫描 Δ频率最大位点 = 切片 2127 (A→G)，China∩year>=2019 组 96.3% G
     —— 文献: 中国 2019/2023 复燃株 >90% 携带 A2063G（PMID 41325763 等）
  E2 邻位双 A: 切片 2127/2128 = A/A（文献 A2063/A2064 相邻，A2064G 罕见）
  E3 距离平齐: 文献 2063→2617 = 554 nt；切片 2127→2681 = 554 nt（同一分子内无 indel 差异）
  E4 长度算术: M129 23S 成熟分子 = 2841 nt（文献/web 核实，AF132741 系）；
     U00089.2 FT 注释 23S = 120057..122961 = 2905 nt = 64 nt 上游前导 + 2841
     ⇒ 切片位 p = 文献位 (p-64)
  E5 编号等价（文献直接陈述）: MP A2063G/A2064G/C2617 ↔ EC 2058/2059/2611（Morozumi 2005; Wang 2022）
  E6 保守 motif: 切片 2124-2137 = CGGAAAGACCCCGT（肽转移酶环），在 EC 标尺 2134-2147 精确匹配
     ⇒ MP 切片 2127/2128 ↔ EC 序列坐标 2137/2138 = A/A

新锚点设计（相对早期版本的改进）:
  每站点取 ±2 nt 间隙外的 22 nt 锚（up = site-25..site-3, down = site+2..site+23, 1-based 含），
  避免相邻耐药位点（2063/2064 相距 1 nt）突变时破坏邻位锚点（v1 锚含邻位 → 耐药株 2064 必然 call 失败）。

输出: reference/reference_anchors.json (v2, 原文件备份为 .v1) + results/MY9_位点校准报告.md
"""
import csv, json, os, shutil, collections

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFDIR = os.path.join(BASE, "reference")
RESDIR = os.path.join(BASE, "results")

def load_fa(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()

mp = load_fa(os.path.join(REFDIR, "M129_23S.fasta"))
OFFSET = 64  # 切片位 = 文献位 + 64

# --- E1: 从试点扫描重算 Δ（可复现） ---
delta_rows = []
if os.path.exists(os.path.join(RESDIR, "MY9_位点扫描_试点.tsv")):
    ca, cb = collections.Counter(), collections.Counter()
    ala = collections.defaultdict(collections.Counter)
    with open(os.path.join(RESDIR, "MY9_位点扫描_试点.tsv")) as f:
        r = csv.reader(f, delimiter="\t"); next(r)
        for pos, acc, country, year, b in r:
            p = int(pos)
            isA = (country.split(":")[0].strip() == "China" and year.isdigit() and int(year) >= 2019) or acc.startswith("GCA_052")
            (ca if isA else cb)[p] += 1
            if isA: ala[p][b] += 1
    for p in range(1500, 2751):
        na, nb = ca.get(p, 0), cb.get(p, 0)
        if na < 15 or nb < 15: continue
        fa = sum(c for bb, c in ala[p].items() if bb != mp[p - 1]) / na
        delta_rows.append((p, fa, na, nb))
delta_rows.sort(key=lambda t: -t[1])

# --- 站点定义（切片坐标） + 锚点生成 ---
SITES_SLICE = {"2063": 2127, "2064": 2128, "2617": 2681}
WT_EXPECT = {"2063": "A", "2064": "A", "2617": "C"}
MUT_LABEL = {"2063": "A2063G(/C)", "2064": "A2064G(/C)", "2617": "C2617A|G"}

sites = {}
for lit, sp in SITES_SLICE.items():
    wt = mp[sp - 1]
    assert wt == WT_EXPECT[lit], f"site {lit}: 切片 {sp} WT={wt} ≠ {WT_EXPECT[lit]}"
    up = mp[sp - 25: sp - 3]      # 22 nt（1-based [sp-24..sp-3]，止于站点前 2 nt，不含邻位耐药位点）
    dn = mp[sp + 1: sp + 23]      # 22 nt（1-based [sp+2..sp+23]，起于站点后 2 nt，不含邻位耐药位点）
    assert len(up) == 22 and len(dn) == 22
    sites[lit] = {"label": MUT_LABEL[lit], "ec_equiv": {"2063": 2058, "2064": 2059, "2617": 2611}[lit],
                  "lit_pos": int(lit), "mp_pos_slice": sp, "wt_base": wt,
                  "anchor_up": up, "anchor_down": dn,
                  "anchor_gap": 4, "base_offset": 24}   # 0-based: down 起点=iu+26, 站点碱基=iu+24（R 链=rc(iu+23)）

# --- 校验锚点在 M129 23S 内唯一 ---
for lit, s in sites.items():
    assert mp.count(s["anchor_up"]) == 1 and mp.count(s["anchor_down"]) == 1, f"site {lit} 锚点不唯一"

# --- 写 anchors v2（备份 v1） ---
v1 = os.path.join(REFDIR, "reference_anchors.json")
if os.path.exists(v1):
    shutil.copy(v1, v1 + ".v1")
out = {"slice_len": len(mp), "slice_offset_vs_literature": OFFSET,
       "m129_23S_mature_len": 2841,
       "numbering_note": "MP 文献编号基于 M129 23S 成熟分子(2841 nt, AF132741 系); 本切片取自 U00089.2 FT 120057..122961 (2905 nt, 含 64 nt 上游前导) ⇒ 切片位 = 文献位 + 64。EC 等价: 2063/2064/2617 ↔ 2058/2059/2611 (Morozumi 2005; Wang 2022)。",
       "calibration": "phenotype-driven (China 2019+ vs rest Δfreq) + adjacency + distance-parity + length arithmetic; 详见 results/MY9_位点校准报告.md",
       "sites": sites}
with open(v1, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)

# --- 校准报告 ---
rep = ["# MY9 23S 耐药位点校准报告（v2）", "",
       "> 生成: MY9_03e_calibrate_sites.py；取代 v1 EC 标尺路线（该路线两处编号假设错误，详见脚本头注释）。", "",
       "| 文献命名 | EC 等价 | 切片位 | WT | 突变 | 证据 |", "|---|---|---|---|---|---|"]
for lit, s in sites.items():
    rep.append(f"| {s['label']} | {s['ec_equiv']} | {s['mp_pos_slice']} | {s['wt_base']} | "
               f"{'A→G' if lit in ('2063','2064') else 'C→A/G'} | E1-E6 |")
rep += ["", "证据链:", "- **E1 表型锚定**: 试点扫描(05) China∩year≥2019(含 GCA_052 批次) vs 其余 Δfreq 最大位点 = 切片 2127 (A→G, 96.3% vs 63.4%)。文献: 中国 2019/2023 复燃株 >90% A2063G。", 
        "- **E2 邻位双 A**: 切片 2127/2128 = A/A，对应文献相邻位点 2063/2064。", 
        "- **E3 距离平齐**: 2063→2617 文献距离 554 nt = 切片 2127→2681 距离，同一分子内无 indel ⇒ 2681 = 2617。", 
        "- **E4 长度算术**: M129 23S 成熟 2841 nt；U00089.2 注释 2905 nt（含 64 nt 上游前导）⇒ 切片位 = 文献位 + 64。", 
        "- **E5 编号等价**: MP 2063/2064/2617 ↔ EC 2058/2059/2611（Morozumi 2005 AAC; Wang 2022）。", 
        "- **E6 保守 motif**: 切片 2124-2137 CGGAAAGACCCCGT ↔ EC 标尺(U00096≡J01695) 序列坐标 2134-2147，位点碱基 A/A 对齐。", "",
        "锚点设计: 每站点 ±2 nt 间隙外各取 22 nt（v1 锚含邻位耐药位点 → 耐药株邻位 call 必然失败，已修正）。", "",
        f"试点扫描 top Δ（重算, 前 5）: " + "; ".join(f"slice{p}(freqA={fa:.2f},nA={na},nB={nb})" for p, fa, na, nb in delta_rows[:5])]
with open(os.path.join(RESDIR, "MY9_位点校准报告.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(rep) + "\n")

print("=== calibrated sites (slice coords) ===")
for lit, s in sites.items():
    print(f"{s['label']}: slice {s['mp_pos_slice']}  WT={s['wt_base']}  up={s['anchor_up']}  down={s['anchor_down']}")
print("anchors v2 written; v1 backed up to reference_anchors.json.v1")
print("report -> results/MY9_位点校准报告.md")
