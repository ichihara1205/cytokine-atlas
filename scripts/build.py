#!/usr/bin/env python3
"""data/ から site/ を生成する。フェーズ1の範囲。

    python3 scripts/build.py
    python3 scripts/build.py --skip-validate   # 検証を飛ばす（下書き確認用）
    python3 scripts/build.py --serve           # 生成してローカルサーバを立てる

生成するページ（CLAUDE.md §6 のうちフェーズ1分）:
    /                     分子一覧 ＋ 5×5通信行列 ＋ 絞り込み
    /molecule/<ID>/       分子ページ
    /families/            3軸切り替えのファミリーツリー ＋「名前と実体のズレ」一覧
    /family/<id>/         ファミリーページ
    /cell/<id>/           細胞ページ（通信辺と分化辺は別セクション）
    /todo/                status: todo と inferred の一覧

site/ は毎回作り直す。site/ を手で編集しない（CLAUDE.md §3）。

このスクリプトが守っている約束:
  - assays.qpcr.primers と assays.antibody は status: verified 以外を出力しない（§0.2）
  - config.yaml の include_own_data: false のとき own_data は HTML にも JSON にも出さない（§2）
  - system_closure は producers/receivers から毎回計算する。データ側の手書きは受け付けない（§5）
  - 通信辺（分子由来）と分化辺（細胞由来）を同じ図に混ぜない（§4）
"""

from __future__ import annotations

import argparse
import html
import json
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from validate import (  # noqa: E402
    ROOT,
    Dataset,
    load_dataset,
    run as run_validate,
    walk_status_blocks,
)

SITE = ROOT / "site"

LEVEL_LABEL = {3: "主要", 2: "明確", 1: "文脈依存", 0: "なし"}

CLOSURE = {
    "closed":   ("系内で完結", "この系の中に出し手も受け手もいる"),
    "out_only": ("出るだけ",   "系内で産生されるが、受け手が系内にいない"),
    "in_only":  ("受けるだけ", "系内に受け手はいるが、出し手が系内にいない"),
    "external": ("系外",       "産生も受容も系内で確認されていない"),
}

EFFECT_LABEL = {
    "fibrosis": "線維化",
    "inflammation": "炎症",
    "lipid": "脂質",
    "differentiation": "分化",
}

AXIS_LABEL = {"structure": "構造相同性", "receptor": "受容体共有", "naming": "慣用名・通し番号"}

# 細胞の表示は「総称」を主にする。ただし由来（初代 / 株 / iPS由来）は必ず併記する。
# 同じ HSC でも LX-2（株）と qHSC（初代）ではデータの読み方が変わるため、
# 総称だけを出すと別物を同一視させてしまう。
ORIGIN_LABEL = {"primary": "初代", "line": "株", "ipsc_derived": "iPS由来"}
STAGE_LABEL = {"pluripotent": "多能性", "progenitor": "前駆", "mature": "成熟"}
STATE_LABEL = {"activated": "活性化", "quiescent": "静止期"}

CLOSURE_RANK = {"closed": 0, "out_only": 1, "in_only": 2, "external": 3}

# 経路のステップ種別。上流→下流の順に色を変える
STEP_LABEL = {"ligand": "リガンド", "receptor": "受容体", "transducer": "伝達",
              "tf": "転写因子", "output": "出力"}
STEP_FILL = {"ligand": "#e8f1fb", "receptor": "#e7f5ec", "transducer": "#f3eefb",
             "tf": "#fdf0e3", "output": "#fdeaea"}
STEP_EDGE = {"ligand": "#5b8fc9", "receptor": "#4aa373", "transducer": "#8e6fc0",
             "tf": "#c98a3d", "output": "#c96a6a"}

# 分子間辺の語彙。validate.py の RELATION_TYPES と対応させる
REL_LABEL = {
    "induces": "誘導する", "antagonizes": "拮抗する", "sequesters": "捕捉する",
    "heterodimer_with": "二量体を作る", "synergizes_with": "相乗する", "opposes": "逆向きに働く",
}
REL_ARROW = {
    "induces": "→", "antagonizes": "⊣", "sequesters": "⊣",
    "heterodimer_with": "+", "synergizes_with": "&", "opposes": "⇄",
}
# 向きのない関係（逆向きの辺を作らない）
REL_SYMMETRIC = {"heterodimer_with", "synergizes_with", "opposes"}
# 分子ページでの言い回し。順方向と逆方向で助詞が変わるので別に持つ
REL_PHRASE = {
    "induces": "{} を誘導する", "antagonizes": "{} に拮抗する",
    "sequesters": "{} を捕捉する", "heterodimer_with": "{} と二量体を作る",
    "synergizes_with": "{} と相乗する", "opposes": "{} と逆向きに働く",
}
REL_PHRASE_REV = {
    "induces": "{} によって誘導される", "antagonizes": "{} に拮抗される",
    "sequesters": "{} に捕捉される",
}


# --------------------------------------------------------------------------
# 小道具
# --------------------------------------------------------------------------
def e(s) -> str:
    return html.escape(str(s), quote=True) if s is not None else ""


EFF_CLS = {-2: "vm2", -1: "vm1", 0: "v0", 1: "vp1", 2: "vp2"}


def eff_cls(v: int) -> str:
    """効果の値をクラス名に。負=青 / 0=灰 / 正=赤 の色分けに使う。"""
    return EFF_CLS.get(v, "v0")


def sign(v: int) -> str:
    """0 は '±0'。明示的に 0 と判断されたことを '+0' より読みやすく示す。"""
    return "±0" if v == 0 else f"{v:+d}"


def status_of(node) -> str | None:
    return node.get("status") if isinstance(node, dict) else None


# --- 細胞の表示名 ---------------------------------------------------------
def cell_of(ds: "Dataset", cid: str) -> dict:
    return ds.cells.get(cid, ({}, ""))[0]


def cell_generic(ds: "Dataset", cid: str) -> str:
    """総称。Hepatocyte / HSC / LSEC / KC …"""
    c = cell_of(ds, cid)
    return c.get("label_generic") or c.get("label_short") or c.get("label_ja") or cid


def cell_detail(ds: "Dataset", cid: str) -> str:
    """由来・実体・状態。'株 / LX-2 / 活性化' のような併記用の短い文字列。"""
    c = cell_of(ds, cid)
    bits = [x for x in (
        ORIGIN_LABEL.get(c.get("origin") or ""),
        c.get("label_short"),
        STATE_LABEL.get(c.get("state") or ""),
    ) if x]
    return " / ".join(bits) if bits else "参照ノード（未培養）"


def cell_color(ds: "Dataset", cid: str) -> str:
    return cell_of(ds, cid).get("color") or "#999"


def cell_dot(ds: "Dataset", cid: str) -> str:
    return (f'<span class="dot" style="background:{e(cell_color(ds, cid))}" '
            f'title="{e(cell_generic(ds, cid))}（{e(cell_detail(ds, cid))}）"></span>')


def render_cell_legend(ctx: "Ctx", ds: "Dataset") -> str:
    """色と細胞の対応表。表の中のドットだけでは何色が何か分からなくなるため、
    一覧のすぐ上に置いて、表をスクロールしても視界に残るようにする。"""
    def chips(ids):
        return "".join(
            f'<a class="legend-chip" href="{ctx.url("cell", cid)}">'
            f'{cell_dot(ds, cid)}<b>{e(cell_generic(ds, cid))}</b>'
            f'<span class="legend-origin">{e(cell_detail(ds, cid))}</span></a>'
            for cid in ids)

    inside = [c for c in ds.cells if cell_of(ds, c).get("in_system", True)]
    outside = [c for c in ds.cells if not cell_of(ds, c).get("in_system", True)]
    return f"""
<div class="legend">
  <div class="legend-row"><span class="legend-head">系内</span>{chips(inside)}</div>
  <div class="legend-row"><span class="legend-head out">系外</span>{chips(outside)}
    <span class="legend-note">系外＝この共培養系にいない。受け手がここにしかない分子は
      「出るだけ」になる</span></div>
</div>"""


def badge(st: str | None) -> str:
    if st == "figure_read":
        return ('<span class="badge badge-figure" title="文献の図から読み取った。'
                '本文に明記がないため目視解釈に依存する">図から読取</span>')
    if st == "inferred":
        return '<span class="badge badge-inferred" title="Claude が一般知識から書いた。未検証">未検証</span>'
    if st == "todo":
        return '<span class="badge badge-todo">TODO</span>'
    return ""


def wrap_status(st: str | None, inner: str) -> str:
    cls = {"inferred": "is-inferred", "todo": "is-todo",
           "figure_read": "is-figure"}.get(st or "", "")
    return f'<div class="block {cls}">{inner}{badge(st)}</div>'


def src_link(src) -> str:
    """出典を表示用リンクにする。PMID / DOI / URL。"""
    if not src or not str(src).strip():
        return ""
    out = []
    for one in [s.strip() for s in str(src).replace(";", ",").split(",") if s.strip()]:
        if one.upper().startswith("PMID:"):
            pid = one.split(":", 1)[1].strip()
            out.append(f'<a class="src" href="https://pubmed.ncbi.nlm.nih.gov/{e(pid)}/" '
                       f'target="_blank" rel="noopener">PMID:{e(pid)}</a>')
        elif one.lower().startswith("http"):
            out.append(f'<a class="src" href="{e(one)}" target="_blank" rel="noopener">出典</a>')
        else:
            doi = one[4:].strip() if one.lower().startswith("doi:") else one
            out.append(f'<a class="src" href="https://doi.org/{e(doi)}" '
                       f'target="_blank" rel="noopener">{e(doi)}</a>')
    return " ".join(out)


class Ctx:
    """base_url を意識したリンク生成。"""

    def __init__(self, base: str):
        self.base = base.rstrip("/")

    def url(self, *parts: str) -> str:
        p = "/".join(str(x).strip("/") for x in parts if str(x).strip("/"))
        return f"{self.base}/{p}/" if p else f"{self.base}/"

    def asset(self, name: str) -> str:
        return f"{self.base}/assets/{name}"


# --------------------------------------------------------------------------
# 導出
# --------------------------------------------------------------------------
def in_system_set(ds: Dataset) -> set[str]:
    listed = (ds.config.get("system") or {}).get("in_system")
    if listed:
        return set(listed)
    return {cid for cid, (c, _) in ds.cells.items() if c.get("in_system", True)}


def edges(block, min_level: int = 1) -> dict[str, dict]:
    """producers / receivers を {cell_id: entry} に正規化。level<min は落とす。"""
    if not isinstance(block, dict):
        return {}
    out = {}
    for cid, v in block.items():
        if not isinstance(v, dict):
            continue
        lv = v.get("level")
        if isinstance(lv, int) and not isinstance(lv, bool) and lv >= min_level:
            out[cid] = v
    return out


def derive_closure(m: dict, inside: set[str]) -> dict:
    """CLAUDE.md §5: system_closure は必ずここで計算する。データ側には書かせない。"""
    prod = edges(m.get("producers"))
    recv = edges(m.get("receivers"))
    p_in = {c for c in prod if c in inside}
    r_in = {c for c in recv if c in inside}
    p_out = {c for c in prod if c not in inside}
    r_out = {c for c in recv if c not in inside}
    if p_in and r_in:
        key = "closed"
    elif p_in:
        key = "out_only"
    elif r_in:
        key = "in_only"
    else:
        key = "external"
    return {
        "key": key,
        "label": CLOSURE[key][0],
        "desc": CLOSURE[key][1],
        "producers_in": sorted(p_in),
        "receivers_in": sorted(r_in),
        "producers_out": sorted(p_out),
        "receivers_out": sorted(r_out),
    }


def public_molecule(m: dict, include_own: bool) -> dict:
    """出力してよい形に落とす。ここを通っていない dict をレンダリングしない。"""
    out = json.loads(json.dumps(m, ensure_ascii=False, default=str))
    out.pop("system_closure", None)          # 手書きされていても無視する

    assays = out.get("assays") or {}
    qpcr = assays.get("qpcr") or {}
    # §0.2: プライマーと抗体は verified 以外サイトに出さない
    kept_p = [p for p in (qpcr.get("primers") or []) if status_of(p) == "verified"]
    dropped_p = len(qpcr.get("primers") or []) - len(kept_p)
    kept_a = [a for a in (assays.get("antibody") or []) if status_of(a) == "verified"]
    dropped_a = len(assays.get("antibody") or []) - len(kept_a)
    if qpcr:
        qpcr["primers"] = kept_p
        assays["qpcr"] = qpcr
    if "antibody" in assays:
        assays["antibody"] = kept_a
    out["_dropped"] = {"primers": dropped_p, "antibody": dropped_a}
    if assays:
        out["assays"] = assays

    if not include_own:
        out.pop("own_data", None)            # §2: 隠すのではなく書き出さない
    return out


def work_view(m: dict, include_own: bool) -> dict:
    """作業リスト用の見え方。

    public_molecule と違い、verified でないプライマー・抗体を落とさない。
    落としてしまうと「未検証だから出力しない」ものが /todo からも消え、
    作業リストとして機能しなくなる。
    ただし own_data の除外だけは公開設定に従う。
    /todo が表示するのは path と note だけなので、配列そのものは出力されない。
    """
    out = json.loads(json.dumps(m, ensure_ascii=False, default=str))
    if not include_own:
        out.pop("own_data", None)
    return out


def cell_coverage(ds: Dataset, mols: dict) -> list[dict]:
    """細胞ごとに producers / receivers に何件書かれているかを数える。

    件数が少ない細胞は「通信が少ない細胞」ではなく「まだ書いていない細胞」である可能性が高い。
    空欄と『調べた上でなし』が見た目で区別できないので、ここで明示的に集計して晒す。"""
    rows = []
    n_mol = len(mols)
    for cid in ds.cells:
        # 「評価済み」= その細胞が producers/receivers にキーとして現れる（level 0 を含む）。
        # level 0 は『調べた上でなし』なので、書いていないことと区別して数える。
        assessed = sum(1 for m in mols.values()
                       if cid in (m.get("producers") or {})
                       or cid in (m.get("receivers") or {}))
        p = sum(1 for m in mols.values() if cid in edges(m.get("producers")))
        r = sum(1 for m in mols.values() if cid in edges(m.get("receivers")))
        rows.append({"id": cid, "producers": p, "receivers": r, "total": p + r,
                     "assessed": assessed, "unassessed": n_mol - assessed,
                     "n_mol": n_mol,
                     "in_system": bool(cell_of(ds, cid).get("in_system", True))})
    rows.sort(key=lambda x: (-x["assessed"], -x["total"], x["id"]))
    return rows


def collect_todo(ds: Dataset, include_own: bool) -> list[dict]:
    rows = []
    stores = (("molecule", ds.molecules), ("cell", ds.cells), ("family", ds.families))
    for kind, store in stores:
        for oid, (obj, relfile) in store.items():
            src = obj
            if kind == "molecule":
                src = work_view(obj, include_own)
            blocks: list = []
            walk_status_blocks(src, "", blocks, set())
            for path, wrapped in blocks:
                node = wrapped["_node"]
                st = node.get("status")
                if st in ("todo", "inferred", "figure_read"):
                    label = (node.get("note") or node.get("text") or node.get("name")
                             or node.get("primary") or node.get("defining_feature") or "")
                    rows.append({
                        "kind": kind, "id": oid, "file": relfile,
                        "path": path or "(top)", "status": st,
                        "label": str(label)[:120],
                        "figure": str(node.get("figure") or ""),
                        "source": str(node.get("source") or ""),
                    })
    order = {"todo": 0, "figure_read": 1, "inferred": 2}
    rows.sort(key=lambda r: (order[r["status"]], r["kind"], r["id"], r["path"]))
    return rows


def build_matrix(ds: Dataset, mols: dict, cell_order: list[str], min_level: int) -> dict:
    """行=産生側 / 列=受容側。セルに分子 id を並べる。"""
    grid = {p: {r: [] for r in cell_order} for p in cell_order}
    for mid, m in mols.items():
        prod = edges(m.get("producers"), min_level)
        recv = edges(m.get("receivers"), min_level)
        for p in prod:
            if p not in grid:
                continue
            for r in recv:
                if r in grid[p]:
                    grid[p][r].append(mid)
    for p in grid:
        for r in grid[p]:
            grid[p][r].sort()
    return grid


# --------------------------------------------------------------------------
# ページ骨格
# --------------------------------------------------------------------------
def page(ctx: Ctx, title: str, body: str, cfg: dict, active: str = "") -> str:
    site = cfg.get("site") or {}
    nav = [("", "分子一覧"), ("pathway", "受容経路"), ("production", "産生経路"),
           ("architecture", "受容体の形"), ("relations", "関係"),
           ("families", "ファミリー"), ("todo", "TODO")]
    links = "".join(
        f'<a href="{ctx.url(slug)}" class="{"active" if active == slug else ""}">{e(label)}</a>'
        for slug, label in nav)
    own = "" if cfg.get("include_own_data") else (
        '<span class="flag" title="config.yaml の include_own_data が false。'
        '未発表の実測は出力されていない">own_data 非出力</span>')
    return f"""<!doctype html>
<html lang="{e(site.get('lang', 'ja'))}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="robots" content="noindex, nofollow">
<title>{e(title)} | {e(site.get('title', 'Cytokine Atlas'))}</title>
<link rel="stylesheet" href="{ctx.asset('style.css')}">
</head>
<body>
<header class="topbar">
  <a class="brand" href="{ctx.url()}">{e(site.get('title', 'Cytokine Atlas'))}</a>
  <nav>{links}</nav>
  {own}
</header>
<main>
{body}
</main>
<footer>
  <p>data/ から自動生成。<strong>site/ を直接編集しない。</strong>
     直したいものがあれば data/ か scripts/ を直して再ビルドする。</p>
</footer>
<script src="{ctx.asset('app.js')}"></script>
</body>
</html>
"""


def write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


# --------------------------------------------------------------------------
# 各ページ
# --------------------------------------------------------------------------
def render_index(ctx: Ctx, ds: Dataset, mols: dict, closures: dict, cfg: dict) -> str:
    cell_order = (cfg.get("matrix") or {}).get("cells") or []
    min_level = (cfg.get("matrix") or {}).get("min_level", 1)
    grid = build_matrix(ds, mols, cell_order, min_level)

    def axis_head(cid: str) -> str:
        return (f'{cell_dot(ds, cid)}<a href="{ctx.url("cell", cid)}">'
                f'{e(cell_generic(ds, cid))}</a>'
                f'<span class="origin">{e(cell_detail(ds, cid))}</span>')

    # --- 行列 ---
    head = "".join(f'<th>{axis_head(c)}</th>' for c in cell_order)
    # 濃淡は実際の最大値に合わせる。固定の閾値だと分子が増えた途端に全部同じ色になる。
    peak = max((len(grid[p][r]) for p in cell_order for r in cell_order), default=0)
    rows = []
    for p in cell_order:
        tds = []
        for r in cell_order:
            ids = grid[p][r]
            if not ids:
                tds.append('<td class="m0"></td>')
                continue
            names = ", ".join(mols[i].get("symbol") or i for i in ids)
            heat = max(1, -(-len(ids) * 5 // peak)) if peak else 1
            tds.append(
                f'<td class="mx h{heat}" data-producer="{e(p)}" data-receiver="{e(r)}" '
                f'data-mols="{e(",".join(ids))}" title="{e(names)}" tabindex="0">'
                f'<span class="n">{len(ids)}</span></td>')
        rows.append(f'<tr><th>{axis_head(p)}</th>{"".join(tds)}</tr>')

    matrix_html = f"""
<section class="card">
  <h2>類洞の断面 <span class="sub">行列の細胞が実際にどこにいるか</span></h2>
  {sinusoid_svg(ctx, ds, mols)}
</section>

<section class="card">
  <h2>通信行列 <span class="sub">行 = 産生側 → 列 = 受容側</span></h2>
  <p class="note">セルをクリックすると下の一覧がその組み合わせに絞られる。
     level {min_level} 未満（＝なし・弱い）の辺は載せていない。
     ここに出るのは<strong>別々の細胞どうしの会話</strong>——ある細胞が出した分子を別の細胞が受け取る、
     という同時点の関係（通信辺）。<strong>iPSC→preMac→iKC のような分化の道筋（分化辺）は含まない。</strong>
     矢印の意味が違うので混ぜていない。分化は各細胞ページで見る。</p>
  <div class="scroll-x">
    <table class="matrix">
      <thead><tr><th></th>{head}</tr></thead>
      <tbody>{"".join(rows)}</tbody>
    </table>
  </div>
  <p class="matrix-status" id="matrix-status"></p>
</section>"""

    # --- 絞り込み ---
    fam_opts = {ax: sorted(
        {(fid, ds.families[fid][0].get("label_ja") or fid)
         for m in mols.values()
         if (fid := ((m.get("families") or {}).get(ax))) and fid in ds.families})
        for ax in AXIS_LABEL}
    fam_selects = "".join(
        f'<label>{e(AXIS_LABEL[ax])}<select data-filter="fam-{ax}">'
        f'<option value="">すべて</option>' +
        "".join(f'<option value="{e(fid)}">{e(lab)}</option>' for fid, lab in opts) +
        "</select></label>"
        for ax, opts in fam_opts.items())

    closure_opts = "".join(
        f'<option value="{k}">{e(v[0])}</option>' for k, v in CLOSURE.items())

    filters = f"""
<section class="card">
  <h2>絞り込み</h2>
  <div class="filters">
    <label>検索<span class="sub">id・symbol・gene・別名・ファミリー名</span>
      <input type="search" data-filter="q" placeholder="IL6 / IL-6 / BSF-2 / gp130"></label>
    {fam_selects}
    <label>系内での閉じ方<select data-filter="closure"><option value="">すべて</option>{closure_opts}</select></label>
    <label>検証状態<select data-filter="status">
      <option value="">すべて</option><option value="verified">verified を含む</option>
      <option value="has-inferred">未検証を含む</option><option value="has-todo">TODO を含む</option>
    </select></label>
    <button type="button" id="reset">解除</button>
  </div>
</section>"""

    # --- 一覧（リスト表示・列見出しクリックで並び替え）---
    def side_cell(m, key):
        """産生／受容を、件数＋細胞色のドットで表す。"""
        ids = sorted(edges(m.get(key)))
        dots = "".join(cell_dot(ds, c) for c in ids)
        names = "、".join(f"{cell_generic(ds, c)}（{cell_detail(ds, c)}）" for c in ids)
        return (f'<td class="num" data-v="{len(ids)}" title="{e(names)}">'
                f'<b>{len(ids)}</b><span class="dots">{dots}</span></td>')

    rows_html = []
    for mid, m in sorted(mols.items(), key=lambda kv: (kv[1].get("symbol") or kv[0])):
        cl = closures[mid]
        eff = m.get("effects") or {}
        blocks: list = []
        walk_status_blocks(m, "", blocks, set())
        n_inf = sum(1 for _, w in blocks if w["_node"].get("status") == "inferred")

        eff_tds = ""
        for k in EFFECT_LABEL:
            v = eff.get(k)
            if isinstance(v, int) and not isinstance(v, bool):
                # バッジは span に入れる。td 自体に display:inline-block が付くと
                # そのセルがテーブルレイアウトから外れて列がずれる。
                eff_tds += (f'<td class="num" data-v="{v}">'
                            f'<span class="effv {eff_cls(v)}">{sign(v)}</span></td>')
            else:
                eff_tds += '<td class="num" data-v="-99"><span class="na">–</span></td>'

        rows_html.append(
            f'<tr class="mol" data-id="{e(mid)}">'
            f'<td class="sym"><a href="{ctx.url("molecule", mid)}">'
            f'{e(m.get("symbol") or mid)}</a></td>'
            f'<td class="gene-col">{e(m.get("gene") or "")}</td>'
            f'<td data-v="{CLOSURE_RANK[cl["key"]]}">'
            f'<span class="closure c-{e(cl["key"])}" title="{e(cl["desc"])}">'
            f'{e(cl["label"])}</span></td>'
            + side_cell(m, "producers") + side_cell(m, "receivers")
            + eff_tds
            + f'<td class="num" data-v="{n_inf}">{n_inf}</td></tr>')

    cols = [("分子", "sym", "text"), ("gene", "gene", "text"),
            ("閉じ方", "closure", "num"), ("産生", "prod", "num"), ("受容", "recv", "num")]
    cols += [(EFFECT_LABEL[k], k, "num") for k in EFFECT_LABEL]
    cols += [("未検証", "inferred", "num")]
    head_html = "".join(
        f'<th class="sortable{" num" if t == "num" else ""}" data-col="{i}" '
        f'data-type="{t}" tabindex="0">{e(lab)}<span class="arrow"></span></th>'
        for i, (lab, _key, t) in enumerate(cols))

    empty = "" if mols else """
  <p class="empty">分子がまだ1件もない。<code>data/molecules/_TEMPLATE.yaml</code> をコピーして
     <code>data/molecules/IL6.yaml</code> のように作ると、ここに並ぶ。</p>"""

    return f"""
<h1>分子一覧</h1>
<p class="lede">産生↔受容マップ。<strong>受け手が系内にいない分子</strong>を見つけることがこのページの主目的。</p>
<div class="disclaimer">
  <b>このサイトの記述の大半は未検証です。</b>
  一次文献または実測で確認できたものだけ <span class="badge badge-inferred"
  style="background:#e7f5ec;color:#166534;border-style:solid">verified</span>
  とし、それ以外は<span class="badge badge-inferred">未検証</span>を付けています。
  未検証の記述は一般知識から書いたもので、裏を取っていません。
  <span class="badge badge-figure">図から読取</span> は文献の図を目視で読み取ったもので、
  本文に明記された verified より根拠が弱いものです。
  内訳と残作業は <a href="{ctx.url("todo")}">TODO</a> にあります。
  数値・配列・品番を引くときは必ず出典に当たってください。
</div>
{matrix_html}
{filters}
<section class="card">
  <h2>分子 <span class="sub" id="count">{len(mols)}</span>
      <span class="sub">　列見出しをクリックで並び替え</span></h2>
  {empty}
  <div class="mol-scroll">
    {render_cell_legend(ctx, ds)}
    <table class="mol-table" id="mol-table">
      <thead><tr>{head_html}</tr></thead>
      <tbody id="mol-grid">{"".join(rows_html)}</tbody>
    </table>
  </div>
  <p class="note">産生・受容の数字は level 1 以上の細胞数。ドットの色は上の対応表の細胞。
     効果は -2 〜 +2、「–」は未記載。</p>
</section>"""


def render_molecule(ctx: Ctx, ds: Dataset, mid: str, m: dict, cl: dict, cfg: dict) -> str:
    def cell_link(cid):
        c = ds.cells.get(cid, ({}, ""))[0]
        tag = "" if c.get("in_system", True) else '<span class="outside">系外</span>'
        return (f'<a class="cellchip" href="{ctx.url("cell", cid)}">'
                f'{cell_dot(ds, cid)}{e(cell_generic(ds, cid))}{tag}'
                f'<span class="origin">{e(cell_detail(ds, cid))}</span></a>')

    def side_table(key, title):
        rows = []
        for cid, v in (m.get(key) or {}).items():
            lv = v.get("level")
            rows.append(
                f'<tr class="{"is-inferred" if v.get("status") == "inferred" else ""}">'
                f'<td>{cell_link(cid)}</td>'
                f'<td><span class="lv lv{lv}">{e(lv)}</span> {e(LEVEL_LABEL.get(lv, ""))}</td>'
                f'<td>{src_link(v.get("source"))} {badge(v.get("status"))}</td>'
                f'<td>{e(v.get("note") or "")}</td></tr>')
        if not rows:
            return f"<h3>{e(title)}</h3><p class='empty'>未記載</p>"
        return (f"<h3>{e(title)}</h3><table class='rel'><thead><tr>"
                f"<th>細胞</th><th>level</th><th>出典</th><th>備考</th></tr></thead>"
                f"<tbody>{''.join(rows)}</tbody></table>")

    fams = m.get("families") or {}
    fam_rows = "".join(
        f'<tr><td>{e(AXIS_LABEL[ax])}</td><td>'
        + (f'<a href="{ctx.url("family", fid)}">'
           f'{e(ds.families[fid][0].get("label_ja") or fid)}</a>'
           if (fid := fams.get(ax)) and fid in ds.families else '<span class="empty">未設定</span>')
        + "</td></tr>"
        for ax in AXIS_LABEL)

    rec = "".join(
        wrap_status(status_of(r), f'<b>{e(r.get("name"))}</b>'
                    + (f' <span class="sub">{e(" + ".join(r.get("subunits") or []))}</span>'
                       if r.get("subunits") else "")
                    + (f'<p>{e(r.get("note"))}</p>' if r.get("note") else "")
                    + f' {src_link(r.get("source"))}')
        for r in (m.get("receptor") or []) if isinstance(r, dict))

    dn = m.get("downstream") or {}
    downstream = wrap_status(status_of(dn),
                             f'<b>{e(dn.get("primary") or "")}</b>'
                             + (f'<p>副経路: {e(", ".join(dn.get("secondary") or []))}</p>'
                                if dn.get("secondary") else "")
                             + f' {src_link(dn.get("source"))}') if dn else ""

    form, rng = m.get("form") or {}, m.get("range") or {}
    props = ""
    if form:
        props += wrap_status(status_of(form),
                             f'分泌形式: <b>{e(form.get("secretion"))}</b>'
                             + (f'<p>潜在型の活性化: {e(form.get("latency"))}</p>'
                                if form.get("latency") else "")
                             + f' {src_link(form.get("source"))}')
    if rng:
        props += wrap_status(status_of(rng),
                             f'到達範囲: <b>{e(rng.get("mode"))}</b>'
                             + (f'<p>{e(rng.get("note"))}</p>' if rng.get("note") else "")
                             + f' {src_link(rng.get("source"))}')

    eff = m.get("effects") or {}
    eff_rows = "".join(
        f'<tr><td>{e(EFFECT_LABEL[k])}</td><td class="effcell">'
        f'<span class="bar {eff_cls(v)}"></span>'
        f'<b class="effv {eff_cls(v)}">{sign(v)}</b></td></tr>'
        for k in EFFECT_LABEL
        if isinstance((v := eff.get(k)), int) and not isinstance(v, bool))
    eff_html = ""
    if eff_rows:
        eff_html = wrap_status(
            status_of(eff),
            f'<table class="eff-table"><tbody>{eff_rows}</tbody></table>'
            + (f'<p class="note">{e(eff.get("note"))}</p>' if eff.get("note") else "")
            + f' {src_link(eff.get("source"))}')

    # --- 測定法 ---
    assays = m.get("assays") or {}
    dropped = m.get("_dropped") or {}
    a_parts = []
    for key, label in (("elisa", "ELISA"), ("array", "アレイ")):
        items = assays.get(key) or []
        if items:
            a_parts.append(f"<h3>{label}</h3>" + "".join(
                wrap_status(status_of(i),
                            f'{e(i.get("vendor") or "")} <code>{e(i.get("cat") or "")}</code> '
                            f'{e(i.get("note") or "")} {src_link(i.get("source"))}')
                for i in items if isinstance(i, dict)))
    primers = ((assays.get("qpcr") or {}).get("primers")) or []
    p_rows = "".join(
        f'<tr><td><code>{e(p.get("fwd"))}</code></td><td><code>{e(p.get("rev"))}</code></td>'
        f'<td>{e(p.get("amplicon") or "")}</td><td>{src_link(p.get("source"))}</td></tr>'
        for p in primers)
    a_parts.append("<h3>qPCR プライマー</h3>" + (
        f'<table class="rel"><thead><tr><th>Fwd</th><th>Rev</th><th>amplicon</th><th>出典</th>'
        f'</tr></thead><tbody>{p_rows}</tbody></table>' if p_rows else
        '<p class="empty">verified なプライマーは登録されていない。'
        '<strong>推測で埋めない。</strong>一次情報が取れたときだけ書く。</p>'))
    if dropped.get("primers"):
        a_parts.append(f'<p class="note">※ verified でないプライマー {dropped["primers"]} 件は'
                       f'出力から除外した。</p>')
    abs_ = assays.get("antibody") or []
    a_parts.append("<h3>抗体</h3>" + ("".join(
        wrap_status(status_of(a),
                    f'{e(a.get("target") or "")} clone <code>{e(a.get("clone") or "")}</code> '
                    f'{e(a.get("vendor") or "")} <code>{e(a.get("cat") or "")}</code> '
                    f'{src_link(a.get("source"))}') for a in abs_)
        if abs_ else '<p class="empty">verified な抗体は登録されていない。</p>'))
    if dropped.get("antibody"):
        a_parts.append(f'<p class="note">※ verified でない抗体 {dropped["antibody"]} 件は'
                       f'出力から除外した。</p>')
    rs = assays.get("rnaseq") or {}
    if rs:
        a_parts.append("<h3>RNA-seq</h3>" + wrap_status(
            status_of(rs), f'{e(rs.get("note") or "")} {src_link(rs.get("source"))}'))

    # --- 自分の実測 ---
    own_html = ""
    if cfg.get("include_own_data"):
        od = m.get("own_data") or []
        if od:
            rows = "".join(
                f'<tr><td>{e(o.get("date"))}</td><td>{e(o.get("experiment") or "")}</td>'
                f'<td>{e(o.get("condition") or "")}</td>'
                f'<td class="num">{e(o.get("value"))} {e(o.get("unit") or "")}</td>'
                f'<td>{e(o.get("note") or "")}</td></tr>' for o in od)
            own_html = f"""
<section class="card own">
  <h2>自分の実測 <span class="sub">未発表</span></h2>
  <table class="rel"><thead><tr><th>日付</th><th>実験</th><th>条件</th><th>値</th><th>備考</th>
  </tr></thead><tbody>{rows}</tbody></table>
</section>"""
    # include_own_data: false のときは own_data の存在自体を出さない。
    # 「実測はあるが伏せた」と書くこと自体が未発表データの存在を漏らすため、
    # ヘッダの「own_data 非出力」表示（サイト全体の状態）だけに留める。

    notes = "".join(
        wrap_status(status_of(n), f'{e(n.get("text") or "")} {src_link(n.get("source"))}')
        for n in (m.get("notes") or []) if isinstance(n, dict))
    refs = "".join(
        f'<li>{e(r.get("citation") or "")} '
        + (src_link(f"PMID:{r['pmid']}") if r.get("pmid") else src_link(r.get("doi")))
        + "</li>"
        for r in (m.get("refs") or []) if isinstance(r, dict))

    out_note = ""
    if cl["key"] == "out_only":
        who = "、".join(
            (ds.cells.get(c, ({}, ""))[0].get("label_ja") or c) for c in cl["receivers_out"])
        out_note = (f'<p class="alarm">系内に受け手がいない。'
                    + (f'受け手として登録されているのは系外の {e(who)} のみ。' if who else
                       '受け手が一つも登録されていない。')
                    + '</p>')

    # 経路の中での位置（分子側には書かせず、経路データから導出する）
    pw_rows = "".join(
        f'<tr><td><a href="{ctx.url("pathway", pid)}">'
        f'{e(ds.pathways[pid][0].get("label_ja") or pid)}</a></td>'
        f'<td><span class="pw-tag k-{e(st.get("kind"))}">'
        f'{e(STEP_LABEL.get(st.get("kind"), st.get("kind")))}</span> {e(st.get("label"))}</td>'
        f'<td class="note">{e(ds.pathways[pid][0].get("outcome") or "")[:90]}</td></tr>'
        for pid, st in molecule_pathways(ds).get(mid, []))
    pw_html = (f'<table class="rel"><thead><tr><th>経路</th><th>この分子の位置</th>'
               f'<th>経路の役割</th></tr></thead><tbody>{pw_rows}</tbody></table>'
               if pw_rows else
               '<p class="unwritten">どの経路にもまだ載せていない。'
               '<strong>経路を持たないという意味ではなく、書いていないだけ。</strong></p>')

    # 産生経路（分子側には書かせず data/production/ から導出する）
    prod_id = production_of(ds).get(mid)
    if prod_id:
        pr = ds.production[prod_id][0]
        cav = measurement_caveats(pr)
        prod_html = (
            ('<div class="alarm"><b>測定上の注意</b><ul>'
             + "".join(f'<li>{e(txt)}</li>' for _k, txt in cav) + '</ul></div>'
             if cav else "")
            + f'<p><a href="{ctx.url("production", prod_id)}">'
            f'{e(pr.get("label_ja") or prod_id)}</a></p>'
            + f'<p class="note">{e(str(pr.get("outcome") or "")[:220])}</p>'
            + production_svg(ctx, ds, pr))
    else:
        prod_html = ('<p class="unwritten">産生経路をまだ書いていない。'
                     '<strong>産生の制御が単純だという意味ではない。</strong> '
                     f'<a href="{ctx.url("production")}">産生経路の一覧</a></p>')

    # 分子間辺。逆向きは導出する
    def msym(x):
        return e(ds.molecules[x][0].get("symbol") or x) if x in ds.molecules else e(x)

    def mlink(x):
        return (f'<a href="{ctx.url("molecule", x)}">{msym(x)}</a>'
                if x in ds.molecules else msym(x))
    rel_rows, seen_rel = [], set()
    for ed in relation_edges(ds):
        for me, other, rev in ((ed["src"], ed["dst"], False), (ed["dst"], ed["src"], True)):
            if me != mid or (other, ed["type"], rev) in seen_rel:
                continue
            seen_rel.add((other, ed["type"], rev))
            if rev and ed["type"] in REL_SYMMETRIC:
                continue
            tmpl = (REL_PHRASE_REV if rev else REL_PHRASE).get(
                ed["type"], "{} — " + str(ed["type"]))
            arrow = REL_ARROW.get(ed["type"], "-")
            if rev and ed["type"] not in REL_SYMMETRIC:
                arrow = "←" if arrow == "→" else "⊢"
            rel_rows.append(
                f'<tr><td class="rel-arrow">{e(arrow)}</td>'
                f'<td>{tmpl.format(mlink(other))}</td>'
                f'<td>{e(ed.get("note") or "")}</td>'
                f'<td>{src_link(ed.get("source"))} {badge(ed.get("status"))}</td></tr>')
    shared = sorted(shares_receptor(ds).get(mid, set()))
    shared_html = ("".join(f'{mlink(s)} ' for s in shared)
                   if shared else '<span class="empty">なし</span>')

    return f"""
<h1>{e(m.get("symbol") or mid)} <span class="gene">{e(m.get("gene") or "")}</span></h1>
<p class="lede">{e(" / ".join(m.get("aliases") or []))}</p>
<p class="oli-line">{oligomer_glyph(molecule_oligomer(m))}
   <span>{e(OLIGOMER_LABEL.get(molecule_oligomer(m) or "", "会合状態は未記載"))}</span>
   <span class="sub">{e((m.get("form") or {}).get("oligomer_note") or "")}</span>
   {src_link((m.get("form") or {}).get("oligomer_source"))}</p>
<p class="closure big c-{e(cl["key"])}">{e(cl["label"])} — {e(cl["desc"])}</p>
{out_note}

<section class="card">
  <h2>ファミリー <span class="sub">3軸は互いに無関係。1本の木にまとめない</span></h2>
  <table class="rel"><tbody>{fam_rows}</tbody></table>
</section>

<section class="card">
  <h2>受容体と下流</h2>
  {rec or '<p class="empty">未記載</p>'}
  {downstream}
</section>

<section class="card">
  <h2>性質</h2>
  {props or '<p class="empty">未記載</p>'}
</section>

<section class="card">
  <h2>どう作られ、どう外に出るか <span class="sub">産生経路</span></h2>
  {prod_html}
</section>

<section class="card">
  <h2>受け取った細胞で何が起きるか <span class="sub">受容経路</span></h2>
  {pw_html}
</section>

<section class="card">
  <h2>他の分子との関係 <span class="sub">分子間辺。細胞どうしの会話とは別のグラフ</span></h2>
  {f'<table class="rel"><tbody>{"".join(rel_rows)}</tbody></table>' if rel_rows
    else '<p class="empty">登録なし</p>'}
  <h3>受容体を共有する分子 <span class="sub">受容体軸から導出</span></h3>
  <div class="chips">{shared_html}</div>
</section>

<section class="card">
  <h2>産生と受容
      <span class="sub">どの細胞が出して、どの細胞が受け取るか（通信辺）</span></h2>
  {sinusoid_svg(ctx, ds, {mid: m}, mid)}
  {side_table("producers", "産生する細胞")}
  {side_table("receivers", "受け取る細胞")}
</section>

<section class="card">
  <h2>効果 <span class="sub">-2 〜 +2。文献が割れたら 0（平均は取らない）</span></h2>
  {eff_html or '<p class="empty">未記載</p>'}
</section>

<section class="card">
  <h2>測定法</h2>
  {"".join(a_parts)}
</section>

{own_html}

<section class="card">
  <h2>メモ</h2>
  {notes or '<p class="empty">なし</p>'}
</section>

<section class="card">
  <h2>文献</h2>
  {f"<ul class='refs'>{refs}</ul>" if refs else '<p class="empty">なし</p>'}
</section>"""


def render_cell(ctx: Ctx, ds: Dataset, cid: str, c: dict, mols: dict, inside: set[str]) -> str:
    # level>=1 は「辺あり」、level 0 は「評価したが、なし」。
    # キーが無いものは「未評価」で、ここには一切現れない（それが区別できることが大事）。
    produces, receives = [], []
    zero_p, zero_r = [], []
    for mid, m in sorted(mols.items()):
        for key, hit, zero in (("producers", produces, zero_p),
                               ("receivers", receives, zero_r)):
            v = (m.get(key) or {}).get(cid)
            if isinstance(v, dict) and isinstance(v.get("level"), int):
                (hit if v["level"] >= 1 else zero).append((mid, m, v))

    def zero_note(items):
        """『調べた上でなし』を明示する。未記載と同じ扱いにしない。"""
        if not items:
            return ""
        links = "、".join(
            f'<a href="{ctx.url("molecule", mid)}">{e(m.get("symbol") or mid)}</a>'
            for mid, m, _ in items)
        return (f'<p class="assessed-zero">評価したうえで「なし」と判断したもの'
                f'（{len(items)}件）: {links}</p>')

    def mol_rows(items):
        if not items:
            return ('<p class="unwritten">まだ1件も記載していない。'
                    '<strong>「この細胞は分子を出さない／受け取らない」という意味ではない。</strong>'
                    'どの細胞を producers / receivers に書くかは執筆時の都合で偏るため、'
                    '空欄は未調査とみなすこと。</p>')
        return ("<table class='rel'><thead><tr><th>分子</th><th>level</th><th>出典</th></tr>"
                "</thead><tbody>" + "".join(
                    f'<tr><td><a href="{ctx.url("molecule", mid)}">'
                    f'{e(m.get("symbol") or mid)}</a></td>'
                    f'<td><span class="lv lv{v["level"]}">{v["level"]}</span> '
                    f'{e(LEVEL_LABEL.get(v["level"], ""))}</td>'
                    f'<td>{src_link(v.get("source"))} {badge(v.get("status"))}</td></tr>'
                    for mid, m, v in items) + "</tbody></table>")

    # 分化辺（通信辺とは別グラフ）
    diff = c.get("differentiation") or {}
    parent = diff.get("from")
    children = [k for k, (v, _) in ds.cells.items()
                if ((v.get("differentiation") or {}).get("from")) == cid]

    def cl(x):
        return (f'<a href="{ctx.url("cell", x)}">'
                f'{e(cell_generic(ds, x))}<span class="origin">'
                f'{e(cell_detail(ds, x))}</span></a>')

    diff_html = ""
    if parent or children or diff:
        parts = []
        if parent:
            parts.append(f"<p>分化元: {cl(parent)}</p>")
        if children:
            parts.append("<p>分化先: " + "、".join(cl(x) for x in sorted(children)) + "</p>")
        if diff.get("protocol"):
            parts.append(f'<p>プロトコル: {e(diff["protocol"])}</p>')
        if diff.get("markers_gained"):
            parts.append("<p>獲得マーカー: " + "".join(
                f'<code>{e(x)}</code> ' for x in diff["markers_gained"]) + "</p>")
        if diff:
            parts.append(f'<p>{src_link(diff.get("source"))} {badge(status_of(diff))}</p>')
        diff_html = "".join(parts)

    notes = "".join(
        wrap_status(status_of(n), f'{e(n.get("text") or "")} {src_link(n.get("source"))}')
        for n in (c.get("notes") or []) if isinstance(n, dict))

    inside_tag = ("系内" if cid in inside else "系外")
    return f"""
<h1><span class="dot big" style="background:{e(c.get("color") or "#999")}"></span>
    {e(cell_generic(ds, cid))}
    <span class="gene">{e(cell_detail(ds, cid))}</span></h1>
<p class="lede">{e(c.get("label_ja") or "")}<br>
   由来: <strong>{e(ORIGIN_LABEL.get(c.get("origin") or "") or "参照ノード（実際には培養していない）")}</strong>
   ／ 分化段階: {e(STAGE_LABEL.get(c.get("stage") or "") or "-")}
   ／ 系統: {e(c.get("lineage") or "-")}
   ／ <strong>{inside_tag}</strong></p>

<section class="card">
  <h2>類洞の中での位置 <span class="sub">この細胞がどこにいるか</span></h2>
  {sinusoid_svg(ctx, ds, mols, focus_cell=cid)}
</section>

<section class="card">
  <h2>出す分子
      <span class="sub">この細胞が産生し、他の細胞に届く分子（通信辺）</span></h2>
  {mol_rows(produces)}
  {zero_note(zero_p)}
</section>

<section class="card">
  <h2>受ける分子
      <span class="sub">他の細胞が出し、この細胞が受け取る分子（通信辺）</span></h2>
  {mol_rows(receives)}
  {zero_note(zero_r)}
</section>

<section class="card">
  <h2>分化
      <span class="sub">この細胞が何から変わり、何に変わるか＝時間の流れ（分化辺）</span></h2>
  <p class="note">上の2つが「同時点にいる別の細胞との会話」なのに対し、これは
     「同じ細胞が時間とともに変わる道筋」。矢印の意味が違うので別に扱う。</p>
  {diff_html or '<p class="empty">分化の記載なし</p>'}
</section>

<section class="card">
  <h2>メモ</h2>
  {notes or '<p class="empty">なし</p>'}
</section>"""




# --------------------------------------------------------------------------
# 分子の形（会合状態）
# --------------------------------------------------------------------------
OLIGOMER_LABEL = {
    "monomer": "単量体", "homodimer": "ホモ二量体", "homotrimer": "ホモ三量体",
    "heterodimer": "ヘテロ二量体", "multimer": "多量体",
}
# 会合状態ごとの円の配置（cx, cy, 色の別）。単位円系で持ち、描画時に拡大する
OLIGOMER_UNITS = {
    "monomer":     [(0, 0, 0)],
    "homodimer":   [(-0.62, 0, 0), (0.62, 0, 0)],
    "homotrimer":  [(0, -0.66, 0), (-0.62, 0.5, 0), (0.62, 0.5, 0)],
    "heterodimer": [(-0.62, 0, 0), (0.62, 0, 1)],
    "multimer":    [(-1.15, 0.1, 0), (-0.4, -0.5, 0), (0.4, -0.5, 0),
                    (1.15, 0.1, 0), (0, 0.62, 0)],
}


def oligomer_shapes(oli: str | None, cx: float, cy: float, r: float) -> str:
    """会合状態を円の並びで描く。四角い箱では二量体と三量体の違いが出ない。"""
    units = OLIGOMER_UNITS.get(oli or "monomer", OLIGOMER_UNITS["monomer"])
    out = []
    for ux, uy, alt in units:
        out.append(f'<circle cx="{cx + ux * r * 1.15:.1f}" cy="{cy + uy * r * 1.15:.1f}" '
                   f'r="{r:.1f}" class="oli-u{alt}"/>')
    return "".join(out)


def oligomer_glyph(oli: str | None, r: float = 6.0) -> str:
    """HTML に埋める単体のグリフ。"""
    w, h = 46, 30
    label = OLIGOMER_LABEL.get(oli or "monomer", "単量体")
    return (f'<svg class="oli" viewBox="0 0 {w} {h}" role="img" '
            f'aria-label="{e(label)}" xmlns="http://www.w3.org/2000/svg">'
            + oligomer_shapes(oli, w / 2, h / 2, r) + '</svg>')


def molecule_oligomer(m: dict) -> str | None:
    return ((m.get("form") or {}).get("oligomer")) if isinstance(m.get("form"), dict) else None



# --------------------------------------------------------------------------
# 類洞の断面図（通信辺を解剖の上に描く）
# --------------------------------------------------------------------------
# 位置は肝類洞の実際の配置に対応させる。中央が類洞内腔、その両側に有窓の内皮、
# 内皮と肝細胞索のあいだがディッセ腔で、そこに星細胞がいる。
# 分化段階のノード（iPSC / preMac / 肝芽細胞）はここに置かない——
# この図は通信辺のためのもので、分化辺と混ぜない（CLAUDE.md §4）。
# 位置は肝類洞の実際の配置に対応させる。中央が類洞内腔、その両側に有窓の内皮、
# 内皮と肝細胞索のあいだがディッセ腔で、そこに星細胞がいる。
# 分化段階のノード（iPSC / preMac / 肝芽細胞）はここに置かない——
# この図は通信辺のためのもので、分化辺と混ぜない（CLAUDE.md §4）。
# 3つ目の値はラベルを置く向き。帯の境界線と衝突するものは右に逃がす。
SINUS_POS = {
    "hepatocyte_pxb": (150, 64, "below"),
    "hsc_lx2":        (300, 146, "right"),
    "hsc_qhsc":       (560, 146, "right"),
    "lsec_tmnk1":     (150, 190, "right"),
    "lsec_ilsec":     (630, 190, "right"),
    "kc_ikc":         (240, 268, "below"),
    "kc_pkc":         (400, 268, "below"),
    "mac_thp1":       (560, 268, "below"),
    "neutrophil":     (90, 268, "below"),
    "t_nk":           (700, 268, "below"),
    "adipocyte":      (612, 478, "right"),
}
SINUS_W, SINUS_H = 760, 512
MOL_HUB = (400, 210)          # 分子そのものを内腔の中央に置く


def sinusoid_svg(ctx: Ctx, ds: Dataset, mols: dict,
                 mid: str | None = None, focus_cell: str | None = None) -> str:
    """肝類洞の断面図。mid を渡すと 産生細胞 → 分子 → 受容細胞 を矢印で重ねる。

    産生と受容を総当たりで結ぶと辺が積になって毛玉になるので、分子を中央のハブに置く。
    辺の数が (産生数 × 受容数) から (産生数 + 受容数) に落ちて読めるようになる。
    """
    m = mols.get(mid) if mid else None
    prod = set(edges(m.get("producers"))) if m else set()
    recv = set(edges(m.get("receivers"))) if m else set()

    g: list[str] = []
    g.append(f'<rect x="0" y="0" width="{SINUS_W}" height="{SINUS_H}" class="sn-bg"/>')
    for y0, y1, cls in ((14, 118, "sn-hep"), (118, 170, "sn-disse"),
                        (170, 322, "sn-lumen"), (322, 374, "sn-disse"),
                        (374, 452, "sn-hep"), (452, 508, "sn-elsewhere")):
        g.append(f'<rect x="0" y="{y0}" width="{SINUS_W}" height="{y1 - y0}" class="{cls}"/>')
    for y, hh in ((18, 92), (380, 66)):
        for x in range(12, SINUS_W - 30, 96):
            g.append(f'<rect x="{x}" y="{y}" width="88" height="{hh}" rx="11" class="sn-hepcell"/>')
    for y in (170, 322):        # 有窓の内皮。破線が fenestrae
        g.append(f'<line x1="0" y1="{y}" x2="{SINUS_W}" y2="{y}" class="sn-endo"/>')
    g.append(f'<text x="10" y="36" class="sn-zone">肝細胞索</text>'
             f'<text x="10" y="136" class="sn-zone">ディッセ腔</text>'
             f'<text x="10" y="478" class="sn-zone">系外の組織（血流の先）</text>'
             f'<text x="10" y="192" class="sn-zone">類洞内腔</text>'
             f'<text x="10" y="206" class="sn-zone2">有窓内皮に囲まれる</text>')

    # --- 通信辺: 産生細胞 → 分子 → 受容細胞 ---
    hx, hy = MOL_HUB
    if m:
        for cid in sorted(prod):
            if cid not in SINUS_POS:
                continue
            x, y, _ = SINUS_POS[cid]
            lv = (m.get("producers") or {}).get(cid, {}).get("level", 1)
            g.append(f'<path d="M{x},{y} Q{(x + hx) / 2:.0f},{(y + hy) / 2 - 26:.0f} '
                     f'{hx - 26},{hy}" class="sn-edge out lv{lv}" marker-end="url(#sn-ar)"/>')
        for cid in sorted(recv):
            if cid not in SINUS_POS:
                continue
            x, y, _ = SINUS_POS[cid]
            lv = (m.get("receivers") or {}).get(cid, {}).get("level", 1)
            g.append(f'<path d="M{hx + 26},{hy} Q{(x + hx) / 2:.0f},{(y + hy) / 2 + 26:.0f} '
                     f'{x},{y}" class="sn-edge in lv{lv}" marker-end="url(#sn-ar)"/>')

    # --- 細胞 ---
    for cid, (cx, cy, side) in SINUS_POS.items():
        if cid not in ds.cells:
            continue
        c = cell_of(ds, cid)
        role = ("both" if cid in prod and cid in recv else
                "prod" if cid in prod else "recv" if cid in recv else "off")
        cls = f'sn-cell r-{role}' + (" is-focus" if cid == focus_cell else "")
        cls += "" if c.get("in_system", True) else " sn-outside"
        lineage = c.get("lineage")
        if lineage == "stellate":
            g.append(f'<path d="M{cx},{cy - 16} l7,10 12,-2 -6,11 6,11 -12,-2 -7,10 '
                     f'-7,-10 -12,2 6,-11 -6,-11 12,2 z" class="{cls}"/>')
        elif lineage == "myeloid":
            g.append(f'<path d="M{cx - 16},{cy} q2,-14 14,-14 q6,-6 15,1 q12,1 11,13 '
                     f'q3,12 -11,13 q-9,6 -16,-1 q-13,0 -13,-12 z" class="{cls}"/>')
        elif lineage == "endothelial":
            g.append(f'<ellipse cx="{cx}" cy="{cy}" rx="23" ry="11" class="{cls}"/>')
        else:
            g.append(f'<rect x="{cx - 22}" y="{cy - 16}" width="44" height="32" rx="8" '
                     f'class="{cls}"/>')
        if side == "right":
            tx, ty, anc = cx + 28, cy - 1, "start"
        else:
            tx, ty, anc = cx, cy + 32, "middle"
        g.append(f'<a href="{ctx.url("cell", cid)}">'
                 f'<text x="{tx}" y="{ty}" class="sn-name" text-anchor="{anc}">'
                 f'{e(cell_generic(ds, cid))}</text>'
                 f'<text x="{tx}" y="{ty + 13}" class="sn-sub" text-anchor="{anc}">'
                 f'{e(c.get("label_short") or "")}</text></a>')

    # --- 中央の分子（会合状態の形で描く）---
    if m:
        oli = molecule_oligomer(m)
        g.append(f'<circle cx="{hx}" cy="{hy}" r="30" class="sn-hub"/>')
        g.append(oligomer_shapes(oli, hx, hy, 9))
        g.append(f'<text x="{hx}" y="{hy + 46}" class="sn-hubname" text-anchor="middle">'
                 f'{e(m.get("symbol") or mid)}</text>'
                 f'<text x="{hx}" y="{hy + 59}" class="sn-sub" text-anchor="middle">'
                 f'{e(OLIGOMER_LABEL.get(oli or "", "会合状態は未記載"))}</text>')

    claim = (f'{m.get("symbol") or mid} を出す細胞と受け取る細胞を肝類洞の断面に重ねた図'
             if m else '肝類洞の断面と、この共培養系の細胞の位置関係')
    return (f'<figure class="sn-fig">'
            f'<svg viewBox="0 0 {SINUS_W} {SINUS_H}" class="sn-svg" role="img" '
            f'aria-label="{e(claim)}" xmlns="http://www.w3.org/2000/svg">'
            '<defs><marker id="sn-ar" viewBox="0 0 10 10" refX="9" refY="5" '
            'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
            '<path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>'
            + "".join(g) + '</svg>'
            f'<figcaption>{e(claim)}。'
            + ("青＝産生側から分子へ、赤＝分子から受容側へ。線の太さが level。" if m else "")
            + '破線の細胞は系外。分化段階のノード（iPSC・前駆マクロファージ・肝芽細胞）は'
            f'時間軸のグラフに属するのでこの図には含めない。</figcaption></figure>')


# --------------------------------------------------------------------------
# 経路（フェーズ2）
# --------------------------------------------------------------------------
def molecule_pathways(ds: Dataset) -> dict[str, list[tuple[str, dict]]]:
    """分子 id -> [(pathway_id, step)]。分子側に pathways を書かせず、ここで導出する。"""
    out: dict[str, list[tuple[str, dict]]] = {}
    for pid, (pw, _) in ds.pathways.items():
        for st in pw.get("steps") or []:
            if not isinstance(st, dict):
                continue
            for sym in st.get("molecules") or []:
                out.setdefault(sym, []).append((pid, st))
    return out


def relation_edges(ds: Dataset) -> list[dict]:
    """分子間辺を1本のリストに集める。

    出どころは3つ:
      1. 分子の relations（手書き）
      2. 経路の output ステップの induces（経路の出力として作られる分子）
      3. 受容体軸のファミリー所属からの導出（shares_receptor_with）— 手で書かせない
    逆向きの辺もここで生成する。
    """
    edges: list[dict] = []
    for mid, (m, _) in ds.molecules.items():
        for rel in m.get("relations") or []:
            if isinstance(rel, dict) and rel.get("target"):
                edges.append({"src": mid, "dst": rel["target"], "type": rel.get("type"),
                              "note": rel.get("note"), "source": rel.get("source"),
                              "status": rel.get("status"), "via": None})
    for pid, (pw, _) in ds.pathways.items():
        for st in pw.get("steps") or []:
            if not isinstance(st, dict) or not st.get("induces"):
                continue
            # その経路のリガンドが、出力として誘導される分子の上流にあたる。
            # ただし枝を跨がせない（TGF-β1 → ヘプシジン のような枝違いの辺ができる）。
            # 拮抗剤のステップも除く（IL-1Ra が IL-6 を誘導することになってしまう）。
            br = st.get("branch")
            ligs = [s for s in (pw.get("steps") or [])
                    if isinstance(s, dict) and s.get("kind") == "ligand"
                    and not s.get("antagonist")
                    and (not br or not s.get("branch") or s.get("branch") == br)]
            srcs = [x for s in ligs for x in (s.get("molecules") or [])]
            for dst in st["induces"]:
                for src in srcs:
                    if src == dst:
                        continue
                    edges.append({"src": src, "dst": dst, "type": "induces",
                                  "note": st.get("note"), "source": pw.get("source"),
                                  "status": pw.get("status"), "via": pid})
    return edges


def shares_receptor(ds: Dataset) -> dict[str, set[str]]:
    """受容体軸のファミリー所属から『受容体を共有する分子』を導出する。
    手で書くと必ず食い違うので、system_closure と同じくここで毎回計算する。"""
    out: dict[str, set[str]] = {}
    by_fam: dict[str, list[str]] = {}
    for mid, (m, _) in ds.molecules.items():
        fid = (m.get("families") or {}).get("receptor")
        if fid:
            by_fam.setdefault(fid, []).append(mid)
    for members in by_fam.values():
        for a in members:
            out.setdefault(a, set()).update(x for x in members if x != a)
    return out


def pathway_svg(ctx: Ctx, ds: Dataset, pw: dict) -> str:
    """steps を細胞の解剖図として描く。手描きせず必ずデータから生成する（CLAUDE.md §2）。

    kind がそのまま細胞内の区画に対応する。箱と矢印の羅列ではなく
    「どこで起きているか」を示すのが目的:
        ligand=細胞外 / receptor=膜を貫通 / transducer=細胞質 / tf=核へ移行 / output=核内
    """
    steps = [s for s in (pw.get("steps") or []) if isinstance(s, dict)]
    by_kind = {k: [s for s in steps if s.get("kind") == k]
               for k in ("ligand", "receptor", "transducer", "tf", "output")}

    W, PAD, BOX_H, ARROW, MEM_H = 760, 18, 58, 34, 16
    NUC_INSET = 34
    body: list[str] = []
    labels: list[str] = []

    def sym(s):
        return (ds.molecules[s][0].get("symbol") or s) if s in ds.molecules else s

    def draw_row(items, y, h, cls, max_bw=330):
        """同じ区画のステップを横並びに置き、中心の x 座標を返す。

        箱を幅いっぱいにしない。受容体が全幅だと膜が隠れて『貫通している』ことが見えず、
        論文図が持つ空間の情報が失われる。
        """
        xs = []
        n = max(len(items), 1)
        inner = W - PAD * 2
        # 1つだけの段は広く取れる。ラベルが切れるより横に伸ばすほうがよい
        bw = min((inner - (n - 1) * 18) / n, 460 if n == 1 else max_bw)
        span = bw * n + 18 * (n - 1)
        x0 = (W - span) / 2
        for i, st in enumerate(items):
            x = x0 + i * (bw + 18)
            body.append(
                f'<rect x="{x:.0f}" y="{y:.0f}" width="{bw:.0f}" height="{h}" rx="8" '
                f'class="pw-box {cls}"/>'
                f'<text x="{x + 11:.0f}" y="{y + 17:.0f}" class="pw-kind">'
                f'{e(STEP_LABEL.get(st.get("kind"), ""))}'
                + (' ・拮抗' if st.get("antagonist") else "") + '</text>'
                f'<text x="{x + 11:.0f}" y="{y + 36:.0f}" class="pw-label">'
                f'{e(str(st.get("label") or "")[:int(bw / 14)])}</text>')
            mols = (st.get("molecules") or []) + (st.get("induces") or [])
            if st.get("kind") == "ligand" and st.get("molecules"):
                olis = {molecule_oligomer(ds.molecules[s][0])
                        for s in st["molecules"] if s in ds.molecules}
                if len(olis) == 1 and (o := olis.pop()):
                    body.append(oligomer_shapes(o, x + bw - 26, y + 20, 5.5))
            if mols:
                cap = max(3, int(bw / 60))
                txt = "、".join(sym(s) for s in mols[:cap])
                if len(mols) > cap:
                    txt += f" ほか{len(mols) - cap}"
                body.append(f'<text x="{x + 11:.0f}" y="{y + 51:.0f}" class="pw-mol">'
                            f'{e(txt[:int(bw / 11)])}</text>')
            xs.append(x + bw / 2)
        return xs

    def connect(top_xs, y0, bot_xs, y1, text):
        """区画をまたぐ矢印。無ラベルの矢印は『なんか関係ある』しか伝えないので必ず注記する。"""
        if not top_xs or not bot_xs:
            return
        pairs = (list(zip(top_xs, bot_xs)) if len(top_xs) == len(bot_xs) > 1
                 else [(a, b) for a in top_xs for b in bot_xs])
        for a, b in pairs:
            body.append(f'<line x1="{a:.0f}" y1="{y0:.0f}" x2="{b:.0f}" y2="{y1:.0f}" '
                        f'class="pw-arrow" marker-end="url(#pw-ar)"/>')
        mid = sum(top_xs) / len(top_xs)
        labels.append(f'<text x="{mid:.0f}" y="{(y0 + y1) / 2 + 4:.0f}" '
                      f'class="pw-edge" text-anchor="middle">{e(text)}</text>')

    # --- 細胞外 ---
    y = 30
    lig_xs = draw_row(by_kind["ligand"], y, BOX_H, "k-ligand")
    y_lig_bottom = y + BOX_H

    # --- 膜（受容体が貫通する）---
    mem_y = y_lig_bottom + ARROW
    rec = by_kind["receptor"]
    rec_y = mem_y - 26
    rec_xs = []
    for i, st in enumerate(rec):
        rcx = W / 2 if len(rec) == 1 else W * (i + 1) / (len(rec) + 1)
        body.append(receptor_glyph(st.get("architecture"), rcx, mem_y, MEM_H))
        arch = arch_label(ds, st.get("architecture"))
        if arch and arch.replace("受容体", "") in str(st.get("label") or ""):
            arch = ""      # ステップ名と同じことを二度書かない
        body.append(f'<text x="{rcx:.0f}" y="{mem_y + MEM_H + 62:.0f}" class="pw-label" '
                    f'text-anchor="middle">{e(str(st.get("label") or "")[:26])}</text>')
        if arch:
            body.append(f'<text x="{rcx:.0f}" y="{mem_y + MEM_H + 76:.0f}" class="pw-mol" '
                        f'text-anchor="middle">{e(arch)}</text>')
        rec_xs.append(rcx)
    connect(lig_xs, y_lig_bottom, rec_xs or [W / 2], rec_y, "結合する")

    y = mem_y + MEM_H + (84 if rec else 26)

    # --- 細胞質 ---
    trans_xs = []
    if by_kind["transducer"]:
        y += ARROW
        trans_xs = draw_row(by_kind["transducer"], y, BOX_H, "k-transducer")
        connect(rec_xs, mem_y + MEM_H + 80, trans_xs, y, "受容体が活性化する")
        y += BOX_H

    # --- 核（tf は境界にまたがせ、output は中に置く）---
    prev_xs, prev_y = (trans_xs or rec_xs), (y if trans_xs else mem_y + MEM_H + 80)
    tf_xs = []
    nuc_top = None
    if by_kind["tf"]:
        y += ARROW
        # 転写因子は核膜の境界にまたがせる。細胞質で活性化してから核へ入るという
        # 順序が、箱の位置そのもので読めるようにするため。
        nuc_top = y + BOX_H / 2
        tf_xs = draw_row(by_kind["tf"], y, BOX_H, "k-tf")
        connect(prev_xs, prev_y, tf_xs, y, "リン酸化して核へ移行")
        prev_xs, prev_y = tf_xs, y + BOX_H
        y += BOX_H

    out_xs = []
    if by_kind["output"]:
        y += ARROW
        # 核は転写因子がある経路にだけ描く。ケモカインの出力は遊走であって転写ではないので、
        # 核を描くと図が事実と食い違う。
        out_xs = draw_row(by_kind["output"], y, BOX_H, "k-output")
        connect(prev_xs, prev_y, out_xs, y,
                "標的遺伝子の転写" if by_kind["tf"] else "細胞の応答")
        y += BOX_H

    H = y + 26
    nuc = ""
    if nuc_top is not None:
        nuc = (f'<rect x="{NUC_INSET}" y="{nuc_top:.0f}" width="{W - NUC_INSET * 2}" '
               f'height="{H - nuc_top - 12:.0f}" rx="26" class="pw-nucleus"/>'
               f'<text x="{NUC_INSET + 12}" y="{nuc_top + 18:.0f}" class="pw-zone">核</text>')

    # 膜（脂質二重層）と区画の地
    mem = [f'<rect x="0" y="0" width="{W}" height="{mem_y:.0f}" class="pw-outside"/>',
           f'<rect x="0" y="{mem_y + MEM_H:.0f}" width="{W}" height="{H - mem_y - MEM_H:.0f}" '
           f'class="pw-cytosol"/>',
           f'<line x1="0" y1="{mem_y:.0f}" x2="{W}" y2="{mem_y:.0f}" class="pw-mem"/>',
           f'<line x1="0" y1="{mem_y + MEM_H:.0f}" x2="{W}" y2="{mem_y + MEM_H:.0f}" '
           f'class="pw-mem"/>']
    for cx in range(8, W, 15):          # 脂質の頭部
        mem.append(f'<circle cx="{cx}" cy="{mem_y + 3:.0f}" r="3" class="pw-lipid"/>'
                   f'<circle cx="{cx}" cy="{mem_y + MEM_H - 3:.0f}" r="3" class="pw-lipid"/>')
    zone = (f'<text x="10" y="18" class="pw-zone">細胞外</text>'
            f'<text x="10" y="{mem_y + MEM_H + 20:.0f}" class="pw-zone">細胞質</text>')

    tail = ("細胞質を経て核で転写を変えるまでの流れ" if by_kind["tf"]
            else "細胞質での応答に至るまでの流れ（この経路は転写ではなく細胞の挙動を変える）")
    claim = f'{pw.get("label_ja") or pw.get("id")}。細胞外のリガンドが膜の受容体に結合し、{tail}'
    return (f'<figure class="pw-fig">'
            f'<svg viewBox="0 0 {W} {H:.0f}" class="pw-svg" role="img" '
            f'aria-label="{e(claim)}" xmlns="http://www.w3.org/2000/svg">'
            '<defs><marker id="pw-ar" viewBox="0 0 10 10" refX="9" refY="5" '
            'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            '<path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>'
            + "".join(mem) + nuc + zone + "".join(body) + "".join(labels)
            + '</svg>'
            f'<figcaption>{e(claim)}</figcaption></figure>')


def arch_label(ds: Dataset, aid: str | None) -> str:
    a = ds.architectures.get(aid or "", ({}, ""))[0]
    return a.get("label_ja") or ""




def receptor_glyph(arch: str | None, cx: float, my: float, mh: float) -> str:
    """受容体を膜での実際の形で描く。

    一律の箱にすると、gp130 が2本鎖であることも、GPCR が膜を7回貫くことも、
    RTK が細胞内にキナーゼドメインを持つことも図から消えてしまう。
    my = 膜の上端 / mh = 膜の厚み。細胞外は上、細胞質は下。
    """
    top, bot = my - 26, my + mh + 26
    g: list[str] = []

    def chain(x, w=13, cls="rc-chain", t=top, b=bot):
        g.append(f'<rect x="{x:.0f}" y="{t:.0f}" width="{w}" height="{b - t:.0f}" '
                 f'rx="5" class="{cls}"/>')

    def kinase(x, y, label=""):
        g.append(f'<rect x="{x:.0f}" y="{y:.0f}" width="26" height="17" rx="4" '
                 f'class="rc-kin"/>')
        if label:
            g.append(f'<text x="{x + 13:.0f}" y="{y + 12:.0f}" class="rc-mini" '
                     f'text-anchor="middle">{e(label)}</text>')

    if arch == "gpcr":
        # 膜を7回貫く。蛇行を7本の柱と上下の連結で表す
        for i in range(7):
            x = cx - 45 + i * 14
            chain(x, 9, "rc-chain", my - 12, my + mh + 12)
            if i < 6:
                y = my - 12 if i % 2 == 0 else my + mh + 12
                g.append(f'<path d="M{x + 4.5:.0f},{y:.0f} q7,{-9 if i % 2 == 0 else 9} '
                         f'14,0" class="rc-loop"/>')
        g.append(f'<circle cx="{cx + 58:.0f}" cy="{my + mh + 22:.0f}" r="10" class="rc-g"/>'
                 f'<text x="{cx + 58:.0f}" y="{my + mh + 26:.0f}" class="rc-mini" '
                 f'text-anchor="middle">Gi</text>')
    elif arch == "tnfr":
        for i in (-1, 0, 1):                       # 三量体
            chain(cx + i * 20 - 6, 12)
    elif arch == "rtk":
        for i in (-1, 1):                          # 二量体化して自己リン酸化
            chain(cx + i * 18 - 6)
            kinase(cx + i * 18 - 13, my + mh + 28, "K")
    elif arch == "stk":
        chain(cx - 26, 13, "rc-chain")             # II型
        chain(cx + 13, 13, "rc-chain2")            # I型
        kinase(cx - 32, my + mh + 28, "II")
        kinase(cx + 7, my + mh + 28, "I")
    elif arch == "tir":
        chain(cx - 20, 13)
        chain(cx + 7, 13, "rc-chain2")
        for i, cy in enumerate((top - 26, top - 12)):   # Ig様ドメイン
            g.append(f'<ellipse cx="{cx - 13:.0f}" cy="{cy:.0f}" rx="11" ry="7" '
                     f'class="rc-ig"/>')
        g.append(f'<rect x="{cx - 26:.0f}" y="{my + mh + 28:.0f}" width="52" height="18" '
                 f'rx="5" class="rc-kin"/>'
                 f'<text x="{cx:.0f}" y="{my + mh + 41:.0f}" class="rc-mini" '
                 f'text-anchor="middle">TIR</text>')
    elif arch == "sefir":
        chain(cx - 7, 14)
        g.append(f'<rect x="{cx - 26:.0f}" y="{my + mh + 28:.0f}" width="52" height="18" '
                 f'rx="5" class="rc-kin"/>'
                 f'<text x="{cx:.0f}" y="{my + mh + 41:.0f}" class="rc-mini" '
                 f'text-anchor="middle">SEFIR</text>')
    else:
        # class1_shared / class2 / 未設定: 共通鎖（濃い）＋固有鎖（淡い）の2本
        chain(cx - 26, 13, "rc-chain2")            # 固有鎖
        chain(cx + 13, 13, "rc-chain")             # 共通鎖
        for i in (-1, 1):                          # 会合する JAK
            g.append(f'<circle cx="{cx + i * 19:.0f}" cy="{my + mh + 36:.0f}" r="9" '
                     f'class="rc-jak"/>'
                     f'<text x="{cx + i * 19:.0f}" y="{my + mh + 40:.0f}" class="rc-mini" '
                     f'text-anchor="middle">JAK</text>')
    return "".join(g)


# --------------------------------------------------------------------------
# 産生経路（中 → 外）
# --------------------------------------------------------------------------
PROD_LABEL = {"trigger": "引き金", "transcription": "転写", "proform": "前駆体",
              "processing": "プロセシング", "release": "放出", "extracellular": "細胞外での活性化"}
PROD_FILL = {"trigger": "#fdf0e3", "transcription": "#f3eefb", "proform": "#eef2f6",
             "processing": "#fdeaea", "release": "#e8f1fb", "extracellular": "#e7f5ec"}
PROD_EDGE = {"trigger": "#c98a3d", "transcription": "#8e6fc0", "proform": "#8b95a1",
             "processing": "#c96a6a", "release": "#5b8fc9", "extracellular": "#4aa373"}
ROUTE_LABEL = {
    "classical": "ER-Golgi を通る古典的分泌", "shedding": "膜型からの切り出し",
    "gasdermin_pore": "ガスダーミンの孔", "cell_death": "細胞死に伴う放出",
    "unconventional": "非古典的経路",
}


def production_of(ds: Dataset) -> dict[str, str]:
    """分子 id -> 産生経路 id。分子側には書かせず導出する。"""
    return {pr.get("molecule"): pid for pid, (pr, _) in ds.production.items()
            if pr.get("molecule")}


def prod_signals(pr: dict) -> list[int]:
    return sorted({st.get("signal") for st in (pr.get("steps") or [])
                   if isinstance(st, dict) and isinstance(st.get("signal"), int)})


def measurement_caveats(pr: dict) -> list[tuple[str, str]]:
    """上清の濃度が転写量を反映しない理由を steps から導出する。

    「2シグナルかどうか」の一つの真偽値では足りない。上清と mRNA が食い違う理由は
    少なくとも3通りあり、原因が違えば対処も違うので分けて出す。
    手で書かず steps から導出する（system_closure と同じ原則）。
    """
    steps = [s for s in (pr.get("steps") or []) if isinstance(s, dict)]
    kinds = [s.get("kind") for s in steps]
    out: list[tuple[str, str]] = []

    # 1. 放出に、転写とは別の引き金がある
    anchor = kinds.index("transcription") if "transcription" in kinds else -1
    if any(k == "trigger" for k in kinds[anchor + 1:]):
        out.append(("release_gated",
                    "放出に転写とは別の引き金がある。転写が上がっても"
                    "その引き金が来なければ外に出ない"))
    # 2. 分泌したあとに活性化が要る
    if "extracellular" in kinds:
        out.append(("post_secretion",
                    "分泌したあと細胞外で活性化される。出た量と効いた量が一致しない"))
    # 3. 膜型が主体で、可溶型は切り出された一部
    if any(s.get("route") == "shedding" for s in steps):
        out.append(("membrane_form",
                    "膜型として働く分があり、上清に出るのは切り出された一部だけ。"
                    "ELISA は膜型の寄与を見ない"))
    # 4. 細胞死に伴って出る
    if any(s.get("route") == "cell_death" for s in steps):
        out.append(("death_release",
                    "放出が細胞死に伴う。上清の値は産生量ではなく死んだ細胞の量を"
                    "反映しうるので、生存率と併せて読む"))
    return out


def has_caveat(pr: dict) -> bool:
    return bool(measurement_caveats(pr))


def production_svg(ctx: Ctx, ds: Dataset, pr: dict) -> str:
    """産生経路を細胞の解剖図で描く。受容側と逆向きに、核から細胞外へ向かって流れる。

    受容側の図は「外から中へ」なので細胞外が上だった。産生側は「中から外へ」なので
    核を上に、膜と細胞外を下に置く。上から下へ読めば、そのまま産生の順序になる。
    """
    steps = [s for s in (pr.get("steps") or []) if isinstance(s, dict)]
    W, PAD, BOX_H, ARROW, MEM_H = 760, 18, 56, 32, 16
    body: list[str] = []
    labels: list[str] = []

    main = [s for s in steps if s.get("kind") != "trigger"]
    triggers = [s for s in steps if s.get("kind") == "trigger"]

    # 本流の縦位置を先に決める
    y = 44
    rows: list[tuple[dict, float]] = []
    mem_y = None
    for st in main:
        if st.get("kind") in ("release", "extracellular") and mem_y is None:
            mem_y = y - 10
            y += MEM_H + 16
        rows.append((st, y))
        y += BOX_H + ARROW
    H = y - ARROW + 30
    if mem_y is None:
        mem_y = H - 40

    # 区画: 上が核、その下が細胞質、下端に膜と細胞外
    body.append(f'<rect x="0" y="0" width="{W}" height="{mem_y:.0f}" class="pd-cytosol"/>')
    body.append(f'<rect x="0" y="{mem_y + MEM_H:.0f}" width="{W}" '
                f'height="{H - mem_y - MEM_H:.0f}" class="pd-outside"/>')
    nuc_rows = [yy for st, yy in rows if st.get("kind") == "transcription"]
    if nuc_rows:
        nx = (W - 380) / 2 - 16
        body.append(f'<rect x="{nx:.0f}" y="{min(nuc_rows) - 22:.0f}" '
                    f'width="{W - nx * 2:.0f}" height="{BOX_H + 42}" rx="26" '
                    f'class="pd-nucleus"/>'
                    f'<text x="{nx + 12:.0f}" y="{min(nuc_rows) - 6:.0f}" '
                    f'class="pd-zone">核</text>')
    body.append(f'<line x1="0" y1="{mem_y:.0f}" x2="{W}" y2="{mem_y:.0f}" class="pw-mem"/>'
                f'<line x1="0" y1="{mem_y + MEM_H:.0f}" x2="{W}" y2="{mem_y + MEM_H:.0f}" '
                f'class="pw-mem"/>')
    for cx in range(8, W, 15):
        body.append(f'<circle cx="{cx}" cy="{mem_y + 3:.0f}" r="3" class="pw-lipid"/>'
                    f'<circle cx="{cx}" cy="{mem_y + MEM_H - 3:.0f}" r="3" class="pw-lipid"/>')
    body.append(f'<text x="10" y="20" class="pd-zone">細胞質</text>'
                f'<text x="10" y="{mem_y + MEM_H + 20:.0f}" class="pd-zone">細胞外</text>')

    # 本流の箱
    BW, x0 = 380, (760 - 380) / 2
    for st, yy in rows:
        k = st.get("kind")
        body.append(
            f'<rect x="{x0:.0f}" y="{yy:.0f}" width="{BW}" height="{BOX_H}" rx="8" '
            f'fill="{PROD_FILL.get(k, "#eee")}" stroke="{PROD_EDGE.get(k, "#999")}" '
            f'stroke-width="1.5"/>'
            f'<text x="{x0 + 12:.0f}" y="{yy + 17:.0f}" class="pw-kind">'
            f'{e(PROD_LABEL.get(k, k))}</text>'
            f'<text x="{x0 + 12:.0f}" y="{yy + 37:.0f}" class="pw-label">'
            f'{e(str(st.get("label") or "")[:26])}</text>')
        extra = ROUTE_LABEL.get(st.get("route") or "", "") or (
            f'プロテアーゼ: {st.get("protease")}' if st.get("protease") else "")
        if extra:
            body.append(f'<text x="{x0 + 12:.0f}" y="{yy + 51:.0f}" class="pw-mol">'
                        f'{e(extra)}</text>')
    for i in range(len(rows) - 1):
        st, yy = rows[i]
        y1 = rows[i + 1][1]
        body.append(f'<line x1="{W / 2:.0f}" y1="{yy + BOX_H:.0f}" x2="{W / 2:.0f}" '
                    f'y2="{y1:.0f}" class="pw-arrow" marker-end="url(#pd-ar)"/>')

    # 引き金は横から刺す。2シグナル制御ならこれが2本になる
    for st in triggers:
        sg = st.get("signal")
        target = None
        for j, (ms, my) in enumerate(rows):
            if sg == 1 and ms.get("kind") == "transcription":
                target = my
            if sg and sg > 1 and ms.get("kind") in ("processing", "release") and target is None:
                target = my
        if target is None:
            target = rows[0][1] if rows else 60
        body.append(
            f'<rect x="14" y="{target + 4:.0f}" width="{x0 - 34:.0f}" height="{BOX_H - 8}" '
            f'rx="8" fill="{PROD_FILL["trigger"]}" stroke="{PROD_EDGE["trigger"]}" '
            f'stroke-width="1.5" stroke-dasharray="5 3"/>'
            f'<text x="24" y="{target + 20:.0f}" class="pw-kind">シグナル {e(sg)}</text>'
            f'<text x="24" y="{target + 38:.0f}" class="pw-label">'
            f'{e(str(st.get("label") or "")[:15])}</text>'
            f'<line x1="{x0 - 18:.0f}" y1="{target + BOX_H / 2:.0f}" x2="{x0 - 2:.0f}" '
            f'y2="{target + BOX_H / 2:.0f}" class="pw-arrow" marker-end="url(#pd-ar)"/>')

    cav = measurement_caveats(pr)
    claim = (f'{pr.get("label_ja")}。核での転写から細胞外へ出るまでの流れ'
             + ("。転写と放出が別々に制御される" if cav else ""))
    return (f'<figure class="pw-fig">'
            f'<svg viewBox="0 0 {W} {H:.0f}" class="pw-svg pd-svg" role="img" '
            f'aria-label="{e(claim)}" xmlns="http://www.w3.org/2000/svg">'
            '<defs><marker id="pd-ar" viewBox="0 0 10 10" refX="9" refY="5" '
            'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            '<path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>'
            + "".join(body) + "".join(labels) + '</svg>'
            f'<figcaption>{e(claim)}。受容側の経路図とは向きが逆で、'
            f'上（核）から下（細胞外）へ読む。</figcaption></figure>')


def render_pathway(ctx: Ctx, ds: Dataset, pid: str, pw: dict, mols: dict) -> str:
    def mol_chip(sym):
        if sym in mols:
            return (f'<a class="chip" href="{ctx.url("molecule", sym)}">'
                    f'{e(mols[sym].get("symbol") or sym)}</a>')
        return f'<span class="chip ghost">{e(sym)}</span>'

    rows = []
    for st in pw.get("steps") or []:
        if not isinstance(st, dict):
            continue
        chips = "".join(mol_chip(s) for s in (st.get("molecules") or []))
        ind = "".join(mol_chip(s) for s in (st.get("induces") or []))
        genes = " ".join(f'<code>{e(g)}</code>' for g in (st.get("genes") or []))
        rows.append(
            f'<tr><td><span class="pw-tag k-{e(st.get("kind"))}">'
            f'{e(STEP_LABEL.get(st.get("kind"), st.get("kind")))}</span></td>'
            f'<td><b>{e(st.get("label"))}</b>'
            + (f'<p class="note">{e(st.get("note"))}</p>' if st.get("note") else "")
            + (f'<div class="chips">{chips}</div>' if chips else "")
            + (f'<div class="chips">出力: {ind}</div>' if ind else "")
            + (f'<p class="sub">{genes}</p>' if genes else "")
            + ((f'<p class="sub">膜での形: <a href="{ctx.url("architecture")}#'
                f'{e(st.get("architecture"))}">{e(arch_label(ds, st.get("architecture")))}</a>'
                f'</p>') if st.get("architecture") else "")
            + "</td></tr>")

    fam = pw.get("family")
    fam_link = ""
    if fam and fam in ds.families:
        fam_link = (f'<p>受容体ファミリー: <a href="{ctx.url("family", fam)}">'
                    f'{e(ds.families[fam][0].get("label_ja") or fam)}</a></p>')
    notes = "".join(
        wrap_status(status_of(n), f'{e(n.get("text") or "")} {src_link(n.get("source"))}')
        for n in (pw.get("notes") or []) if isinstance(n, dict))

    return f"""
<h1>{e(pw.get("label_ja") or pid)}</h1>
<p class="lede">経路図。<a href="{ctx.url("pathway")}">経路の一覧</a></p>

<section class="card">
  <h2>この経路は何をするか</h2>
  {wrap_status(status_of(pw), f'<p class="outcome">{e(pw.get("outcome") or "")}</p>'
               + src_link(pw.get("source")))}
  {fam_link}
</section>

<section class="card">
  <h2>流れ <span class="sub">上から下へ</span></h2>
  <div class="scroll-x">{pathway_svg(ctx, ds, pw)}</div>
</section>

<section class="card">
  <h2>各段階</h2>
  <table class="rel pw-steps"><tbody>{"".join(rows)}</tbody></table>
</section>

<section class="card">
  <h2>メモ</h2>
  {notes or '<p class="empty">なし</p>'}
</section>"""



def render_architectures(ctx: Ctx, ds: Dataset) -> str:
    """受容体の形と、その形が要求する伝達手段を並べて比べる。

    「1回貫通か7回貫通か」「何本鎖か」は見た目の分類ではなく、
    細胞外の結合をどうやって膜の内側へ渡すかという手段を決めてしまう。
    形と機構を並べて置くことでその対応が読める。
    """
    def panel(a: dict, x: float) -> str:
        """膜と受容体を1枚ぶん描く。"""
        my, mh, w = 74, 16, 200
        cx = x + w / 2
        g = [f'<rect x="{x}" y="18" width="{w}" height="{my - 18}" class="pw-outside"/>',
             f'<rect x="{x}" y="{my + mh}" width="{w}" height="{178 - my - mh}" '
             f'class="pw-cytosol"/>',
             f'<line x1="{x}" y1="{my}" x2="{x + w}" y2="{my}" class="pw-mem"/>',
             f'<line x1="{x}" y1="{my + mh}" x2="{x + w}" y2="{my + mh}" class="pw-mem"/>']
        for cxx in range(int(x) + 7, int(x + w), 15):
            g.append(f'<circle cx="{cxx}" cy="{my + 3}" r="3" class="pw-lipid"/>'
                     f'<circle cx="{cxx}" cy="{my + mh - 3}" r="3" class="pw-lipid"/>')
        g.append(receptor_glyph(a.get("id"), cx, my, mh))
        g.append(f'<text x="{cx:.0f}" y="{194}" class="pw-label" text-anchor="middle">'
                 f'{e(a.get("label_ja") or a.get("id"))[:16]}</text>')
        g.append(f'<text x="{cx:.0f}" y="{210}" class="pw-mol" text-anchor="middle">'
                 f'膜貫通 {e(a.get("tm"))} 回 ／ {e(a.get("chains"))} 本鎖'
                 + ("／ 自前のキナーゼあり" if a.get("intrinsic_kinase") else "") + '</text>')
        return "".join(g)

    # 比べたいのは「1回貫通は寄らないと伝わらない」対「7回貫通は1本でねじれて伝える」
    key = [ds.architectures[k][0] for k in ("class1_shared", "rtk", "gpcr")
           if k in ds.architectures]
    W = 200 * len(key) + 20 * (len(key) - 1)
    cmp_svg = (f'<figure class="pw-fig"><svg viewBox="0 0 {W} 224" class="pw-svg" '
               f'role="img" aria-label="膜貫通回数と鎖数が伝達手段をどう決めるかの比較" '
               f'xmlns="http://www.w3.org/2000/svg">'
               + "".join(panel(a, i * 220) for i, a in enumerate(key))
               + '</svg><figcaption>同じ「受容体」でも膜での形が違えば、'
               '膜を越える手段そのものが変わる。左2つは1回しか貫通しないので'
               '寄らないと伝わらない。右は1本の鎖が7回貫いて樽を作り、'
               'そのねじれ自体が細胞内へ伝わる。</figcaption></figure>')

    rows = []
    for aid, (a, _) in sorted(ds.architectures.items(),
                              key=lambda kv: (kv[1][0].get("tm") == 7, kv[0])):
        users = sorted(pid for pid, (pw, _) in ds.pathways.items()
                       for st in (pw.get("steps") or [])
                       if isinstance(st, dict) and st.get("architecture") == aid)
        links = "、".join(
            f'<a href="{ctx.url("pathway", pid)}">'
            f'{e(ds.pathways[pid][0].get("label_ja") or pid)}</a>' for pid in users)
        rows.append(
            f'<tr id="{e(aid)}"><td><b>{e(a.get("label_ja") or aid)}</b>'
            f'<p class="sub">膜貫通 {e(a.get("tm"))} 回 ／ {e(a.get("chains"))} 本鎖'
            + ("／ 自前のキナーゼあり" if a.get("intrinsic_kinase")
               else "／ キナーゼは会合させて借りる") + '</p></td>'
            f'<td>{e(a.get("crossing") or "")}</td>'
            f'<td>{e(a.get("why_chains") or "")}</td>'
            f'<td>{links or "<span class=empty>未使用</span>"}'
            f'<p>{src_link(a.get("source"))} {badge(status_of(a))}</p></td></tr>')

    return f"""
<h1>受容体の形</h1>
<p class="lede">膜貫通回数と鎖数は見た目の分類ではない。
   <strong>細胞外での結合をどうやって膜の内側へ渡すか</strong>という手段を決めてしまうので、
   形と伝達機構は切り離せない。</p>

<section class="card">
  <h2>1回貫通と7回貫通で何が変わるか</h2>
  {cmp_svg}
</section>

<section class="card">
  <h2>一覧</h2>
  <div class="scroll-x">
  <table class="rel arch-table">
    <thead><tr><th>形</th><th>膜を越える手段</th><th>鎖を分ける意味</th><th>使う経路</th></tr></thead>
    <tbody>{"".join(rows)}</tbody>
  </table>
  </div>
</section>"""


def render_production(ctx: Ctx, ds: Dataset, pid: str, pr: dict, mols: dict) -> str:
    mol = pr.get("molecule")
    rows = []
    for st in pr.get("steps") or []:
        if not isinstance(st, dict):
            continue
        k = st.get("kind")
        extra = []
        if st.get("route"):
            extra.append(ROUTE_LABEL.get(st["route"], st["route"]))
        if st.get("protease"):
            extra.append(f'プロテアーゼ: <code>{e(st["protease"])}</code>')
        if isinstance(st.get("signal"), int):
            extra.append(f'シグナル {st["signal"]}')
        rows.append(
            f'<tr><td><span class="pd-tag p-{e(k)}">{e(PROD_LABEL.get(k, k))}</span></td>'
            f'<td><b>{e(st.get("label"))}</b>'
            + (f'<p class="note">{e(st.get("note"))}</p>' if st.get("note") else "")
            + (f'<p class="sub">{" ／ ".join(extra)}</p>' if extra else "")
            + f' {src_link(st.get("source"))} {badge(status_of(st))}</td></tr>')

    notes = "".join(
        wrap_status(status_of(n), f'{e(n.get("text") or "")} {src_link(n.get("source"))}')
        for n in (pr.get("notes") or []) if isinstance(n, dict))
    cav = measurement_caveats(pr)
    warn = ('<div class="alarm"><b>測定上の注意</b><ul>'
            + "".join(f'<li>{e(txt)}</li>' for _k, txt in cav)
            + '</ul></div>' if cav else "")
    mol_link = (f'<a href="{ctx.url("molecule", mol)}">'
                f'{e(mols[mol].get("symbol") or mol)}</a>' if mol in mols else e(mol))
    return f"""
<h1>{e(pr.get("label_ja") or pid)}</h1>
<p class="lede">産生経路。対象分子: {mol_link}　
   <a href="{ctx.url("production")}">産生経路の一覧</a></p>
{warn}

<section class="card">
  <h2>要点 <span class="sub">何が律速か</span></h2>
  {wrap_status(status_of(pr), f'<p class="outcome">{e(pr.get("outcome") or "")}</p>'
               + src_link(pr.get("source")))}
</section>

<section class="card">
  <h2>流れ <span class="sub">核（上）から細胞外（下）へ</span></h2>
  <div class="scroll-x">{production_svg(ctx, ds, pr)}</div>
</section>

<section class="card">
  <h2>各段階</h2>
  <table class="rel pw-steps"><tbody>{"".join(rows)}</tbody></table>
</section>

<section class="card">
  <h2>メモ</h2>
  {notes or '<p class="empty">なし</p>'}
</section>"""


def render_production_index(ctx: Ctx, ds: Dataset, mols: dict) -> str:
    two_step, single = [], []
    for pid, (pr, _) in sorted(ds.production.items()):
        mol = pr.get("molecule")
        sig = prod_signals(pr)
        routes = sorted({ROUTE_LABEL.get(st.get("route"), "")
                         for st in (pr.get("steps") or [])
                         if isinstance(st, dict) and st.get("route")})
        row = (f'<tr><td><a href="{ctx.url("production", pid)}">'
               f'{e(mols[mol].get("symbol") or mol) if mol in mols else e(mol)}</a></td>'
               f'<td>{e(pr.get("label_ja") or pid)}</td>'
               f'<td class="num">{len(sig) or "–"}</td>'
               f'<td>{e("、".join(r for r in routes if r))}</td>'
               f'<td class="note">{e(str(pr.get("outcome") or "")[:110])}</td></tr>')
        (two_step if has_caveat(pr) else single).append(row)

    covered = set(production_of(ds))
    missing = sorted(set(mols) - covered)
    miss_html = "".join(
        f'<a class="chip" href="{ctx.url("molecule", s)}">'
        f'{e(mols[s].get("symbol") or s)}</a>' for s in missing)

    def table(rs):
        return ('<div class="scroll-x"><table class="rel"><thead><tr><th>分子</th>'
                '<th>経路</th><th class="num">シグナル数</th><th>放出の手段</th>'
                '<th>要点</th></tr></thead><tbody>' + "".join(rs) + '</tbody></table></div>'
                if rs else '<p class="empty">なし</p>')

    return f"""
<h1>産生経路</h1>
<p class="lede">サイトカインが<strong>どう作られて、どうやって細胞の外に出るか</strong>。
   受容側の<a href="{ctx.url("pathway")}">経路</a>が「外から中へ」なのに対し、こちらは「中から外へ」。
   矢印の向きが逆なので別のグラフとして扱う。</p>

<section class="card">
  <h2>上清が転写量を反映しない <span class="sub">{len(two_step)} 件</span></h2>
  <p class="note"><strong>この群は qPCR だけで産生を判断できない。</strong>
     理由は「放出に別の引き金が要る」「分泌後に活性化が要る」「膜型が主体」
     「細胞死に伴って出る」の4通りあり、原因が違えば対処も違う。
     各ページに理由を出している。判定は steps から導出（手で書いていない）。</p>
  {table(two_step)}
</section>

<section class="card">
  <h2>制御点が転写側にある <span class="sub">{len(single)} 件</span></h2>
  <p class="note">放出・活性化に独立した制御点を持たない群。
     mRNA と上清が比較的素直に対応する。</p>
  {table(single)}
</section>

<section class="card">
  <h2>産生経路が未記載 <span class="sub">{len(missing)} / {len(mols)}</span></h2>
  <p class="note">まだ書いていない分子。<strong>産生の制御が単純だという意味ではない。</strong></p>
  <div class="chips">{miss_html or '<span class="empty">なし</span>'}</div>
</section>"""


def render_pathway_index(ctx: Ctx, ds: Dataset, mols: dict) -> str:
    inpw = molecule_pathways(ds)
    cards = []
    for pid, (pw, _) in sorted(ds.pathways.items()):
        ligs = [s for s in (pw.get("steps") or [])
                if isinstance(s, dict) and s.get("kind") == "ligand"]
        chips = "".join(
            f'<a class="chip" href="{ctx.url("molecule", s)}">'
            f'{e(mols[s].get("symbol") or s)}</a>'
            for st in ligs for s in (st.get("molecules") or []) if s in mols)
        cards.append(f"""
    <article class="pw-card">
      <h3><a href="{ctx.url("pathway", pid)}">{e(pw.get("label_ja") or pid)}</a>
          {badge(status_of(pw))}</h3>
      <p class="note">{e(str(pw.get("outcome") or "")[:150])}</p>
      <div class="chips">{chips}</div>
    </article>""")

    orphan = sorted(set(mols) - set(inpw))
    orphan_html = "".join(
        f'<a class="chip" href="{ctx.url("molecule", s)}">'
        f'{e(mols[s].get("symbol") or s)}</a>' for s in orphan)
    return f"""
<h1>経路</h1>
<p class="lede">どの分子が経路のどこにいて、その経路が何をするか。
   分子の <code>producers</code>/<code>receivers</code> が細胞どうしの会話（通信辺）なのに対し、
   ここは<strong>分子の内側で何が起きるか</strong>を扱う。</p>
<section class="card">
  <div class="pw-grid">{"".join(cards)}</div>
</section>

<section class="card">
  <h2>経路未割当 <span class="sub">{len(orphan)} / {len(mols)}</span></h2>
  <p class="note">どの経路にもまだ載せていない分子。
     <strong>経路を持たないという意味ではなく、書いていないだけ。</strong></p>
  <div class="chips">{orphan_html or '<span class="empty">なし</span>'}</div>
</section>"""


def render_relations(ctx: Ctx, ds: Dataset, mols: dict) -> str:
    edges = relation_edges(ds)
    shares = shares_receptor(ds)

    def sym(x):
        return e(mols[x].get("symbol") or x) if x in mols else e(x)

    def link(x):
        return f'<a href="{ctx.url("molecule", x)}">{sym(x)}</a>' if x in mols else sym(x)

    seen, rows = set(), []
    for ed in edges:
        key = (ed["src"], ed["dst"], ed["type"])
        if key in seen:
            continue
        seen.add(key)
        via = (f'<a class="src" href="{ctx.url("pathway", ed["via"])}">'
               f'{e(ds.pathways[ed["via"]][0].get("label_ja") or ed["via"])}</a>'
               if ed.get("via") and ed["via"] in ds.pathways else "")
        rows.append(
            f'<tr class="rel-row" data-type="{e(ed["type"])}">'
            f'<td>{link(ed["src"])}</td>'
            f'<td class="rel-arrow">{e(REL_ARROW.get(ed["type"], "-"))}</td>'
            f'<td>{link(ed["dst"])}</td>'
            f'<td><span class="rel-tag">{e(REL_LABEL.get(ed["type"], ed["type"]))}</span></td>'
            f'<td>{e(ed.get("note") or "")}</td>'
            f'<td>{src_link(ed.get("source"))} {badge(ed.get("status"))} {via}</td></tr>')

    opts = "".join(f'<option value="{e(k)}">{e(v)}</option>' for k, v in REL_LABEL.items())
    n_shared = sum(len(v) for v in shares.values()) // 2
    return f"""
<h1>分子どうしの関係</h1>
<p class="lede">これは<strong>3本目のグラフ</strong>。細胞どうしの会話（通信辺）でも、
   細胞の分化（分化辺）でもなく、<strong>分子と分子の間</strong>の関係を集めたもの。</p>

<section class="card">
  <h2>関係 <span class="sub" id="rel-count">{len(rows)}</span></h2>
  <div class="filters">
    <label>種類<select id="rel-filter"><option value="">すべて</option>{opts}</select></label>
  </div>
  <div class="scroll-x">
  <table class="rel">
    <thead><tr><th></th><th></th><th></th><th>関係</th><th>内容</th><th>出典</th></tr></thead>
    <tbody id="rel-body">{"".join(rows) or '<tr><td colspan="6" class="empty">未登録</td></tr>'}</tbody>
  </table>
  </div>
  <p class="note">誘導（→）の多くは<strong>経路の出力から導出したもの</strong>で、出典欄の経路名がその由来。粒度は「その経路のリガンドなら誘導しうる」までであって、リガンド1つずつを個別に確かめた結果ではない。枝を跨ぐ辺と拮抗剤からの辺は除いてある。</p>
</section>

<section class="card">
  <h2>受容体の共有 <span class="sub">{n_shared} 組・導出</span></h2>
  <p class="note">受容体軸のファミリー所属から毎回計算している。
     <strong>手で書かない</strong>——書くとファミリー定義と必ず食い違う（<code>system_closure</code> と同じ理由）。</p>
  <table class="rel"><tbody>{"".join(
      f'<tr><td>{link(a)}</td><td>{"、".join(link(b) for b in sorted(v))}</td></tr>'
      for a, v in sorted(shares.items()) if v)}</tbody></table>
</section>"""


def render_families(ctx: Ctx, ds: Dataset, mols: dict) -> str:
    gene_to_mol = {str(m.get("gene")): mid for mid, m in mols.items() if m.get("gene")}

    def member_chip(sym):
        mid = gene_to_mol.get(sym) or (sym if sym in mols else None)
        if mid:
            return f'<a class="chip" href="{ctx.url("molecule", mid)}">{e(sym)}</a>'
        return f'<span class="chip ghost" title="分子ファイルはまだ無い">{e(sym)}</span>'

    panels = []
    for ax, label in AXIS_LABEL.items():
        fams = {fid: f for fid, (f, _) in ds.families.items() if f.get("axis") == ax}
        roots = [fid for fid, f in fams.items() if not f.get("parent")]

        def node(fid, depth=0):
            f = fams[fid]
            kids = [k for k, v in fams.items() if v.get("parent") == fid]
            body = (f'<div class="fam-node d{depth}">'
                    f'<a class="fam-name" href="{ctx.url("family", fid)}">'
                    f'{e(f.get("label_ja") or fid)}</a>'
                    f'{badge(status_of(f))}'
                    f'<p class="sub">{e(f.get("defining_feature") or "")}</p>'
                    f'<div class="chips">'
                    + "".join(member_chip(s) for s in (f.get("members") or []))
                    + "</div></div>")
            return body + "".join(node(k, depth + 1) for k in sorted(kids))

        panels.append(f'<div class="axis-panel" data-axis="{e(ax)}" hidden>'
                      + ("".join(node(r) for r in sorted(roots))
                         or '<p class="empty">この軸のファミリーはまだ無い</p>')
                      + "</div>")

    tabs = "".join(
        f'<button type="button" class="tab" data-axis="{e(ax)}">{e(lab)}</button>'
        for ax, lab in AXIS_LABEL.items())

    mm_rows = []
    for fid, (f, _) in sorted(ds.families.items()):
        for mm in f.get("mismatch") or []:
            if not isinstance(mm, dict):
                continue
            mm_rows.append(
                f'<tr class="{"is-inferred" if status_of(mm) == "inferred" else ""}">'
                f'<td>{e(AXIS_LABEL.get(f.get("axis"), f.get("axis")))}</td>'
                f'<td><a href="{ctx.url("family", fid)}">{e(f.get("label_ja") or fid)}</a></td>'
                f'<td>{member_chip(mm["member"]) if mm.get("member") else "<span class=sub>全体</span>"}</td>'
                f'<td>{e(mm.get("note") or "")}</td>'
                f'<td>{src_link(mm.get("source"))} {badge(status_of(mm))}</td></tr>')

    return f"""
<h1>ファミリー</h1>
<p class="lede">サイトカインの命名には<strong>互いに無関係な3つの論理</strong>が混在している。
   1本の木では表せないので、3軸を別々に持って切り替える。</p>

<section class="card">
  <div class="tabs">{tabs}</div>
  {"".join(panels)}
</section>

<section class="card">
  <h2>名前と実体のズレ <span class="sub">このレイヤの存在意義</span></h2>
  <p class="note">名前から一族を推測すると間違える箇所を、ここに集めている。</p>
  <div class="scroll-x">
  <table class="rel mismatch">
    <thead><tr><th>軸</th><th>ファミリー</th><th>分子</th><th>ズレの内容</th><th>出典</th></tr></thead>
    <tbody>{"".join(mm_rows) or '<tr><td colspan="5" class="empty">未登録</td></tr>'}</tbody>
  </table>
  </div>
</section>"""


def render_family(ctx: Ctx, ds: Dataset, fid: str, f: dict, mols: dict) -> str:
    gene_to_mol = {str(m.get("gene")): mid for mid, m in mols.items() if m.get("gene")}

    def member_chip(sym):
        mid = gene_to_mol.get(sym) or (sym if sym in mols else None)
        return (f'<a class="chip" href="{ctx.url("molecule", mid)}">{e(sym)}</a>' if mid
                else f'<span class="chip ghost" title="分子ファイルはまだ無い">{e(sym)}</span>')

    par = f.get("parent")
    kids = [k for k, (v, _) in ds.families.items() if v.get("parent") == fid]
    rel = []
    for pid, (pw, _) in sorted(ds.pathways.items()):
        if pw.get("family") == fid:
            rel.append(f'<p>経路: <a href="{ctx.url("pathway", pid)}">'
                       f'{e(pw.get("label_ja") or pid)}</a></p>')
    if par and par in ds.families:
        rel.append(f'<p>上位: <a href="{ctx.url("family", par)}">'
                   f'{e(ds.families[par][0].get("label_ja") or par)}</a></p>')
    if kids:
        rel.append("<p>分岐: " + "、".join(
            f'<a href="{ctx.url("family", k)}">'
            f'{e(ds.families[k][0].get("label_ja") or k)}</a>' for k in sorted(kids)) + "</p>")

    mm = "".join(
        wrap_status(status_of(x),
                    (f'<b>{member_chip(x["member"])}</b> ' if x.get("member") else
                     '<b class="sub">ファミリー全体</b> ')
                    + f'{e(x.get("note") or "")} {src_link(x.get("source"))}')
        for x in (f.get("mismatch") or []) if isinstance(x, dict))

    return f"""
<h1>{e(f.get("label_ja") or fid)}</h1>
<p class="lede">軸: <strong>{e(AXIS_LABEL.get(f.get("axis"), f.get("axis")))}</strong>
   　<a href="{ctx.url("families")}">3軸を見る</a></p>

<section class="card">
  <h2>この一族を定義するもの</h2>
  {wrap_status(status_of(f), f'{e(f.get("defining_feature") or "")} {src_link(f.get("source"))}')}
  {"".join(rel)}
</section>

<section class="card">
  <h2>所属分子 <span class="sub">薄い表示は分子ファイル未作成</span></h2>
  <div class="chips">{"".join(member_chip(s) for s in (f.get("members") or [])) or '<span class="empty">未登録</span>'}</div>
</section>

<section class="card">
  <h2>名前と実体のズレ・例外</h2>
  {mm or '<p class="empty">登録なし</p>'}
</section>"""


def render_todo(ctx: Ctx, ds: Dataset, rows: list[dict], counts: dict,
                coverage: list[dict]) -> str:
    def table(rs):
        if not rs:
            return '<p class="empty">なし</p>'
        return ("<div class='scroll-x'><table class='rel todo'><thead><tr>"
                "<th>種別</th><th>対象</th><th>位置</th><th>内容</th><th>ファイル</th>"
                "</tr></thead><tbody>" + "".join(
                    f'<tr><td>{e(r["kind"])}</td>'
                    f'<td><a href="{ctx.url(r["kind"], r["id"])}">{e(r["id"])}</a></td>'
                    f'<td><code>{e(r["path"])}</code></td>'
                    f'<td>{e(r["label"])}</td>'
                    f'<td><code>{e(r["file"])}</code></td></tr>' for r in rs)
                + "</tbody></table></div>")

    n_mol = coverage[0]["n_mol"] if coverage else 0
    cov_rows = "".join(
        f'<tr class="{"cov-thin" if c["assessed"] < n_mol * 0.5 else ""}">'
        f'<td><a href="{ctx.url("cell", c["id"])}">{e(cell_generic(ds, c["id"]))}</a>'
        f'<span class="origin">{e(cell_detail(ds, c["id"]))}</span></td>'
        f'<td>{"系内" if c["in_system"] else "系外"}</td>'
        f'<td class="num"><b>{c["assessed"]}</b> / {n_mol}</td>'
        f'<td class="num">{c["producers"]}</td><td class="num">{c["receivers"]}</td>'
        f'<td class="num unassessed-n">{c["unassessed"]}</td></tr>'
        for c in coverage)
    cov_html = f"""
<section class="card">
  <h2>細胞ごとの記載状況 <span class="sub">「未評価」と「調べた上でなし」は別物</span></h2>
  <p class="note"><strong>評価済み</strong>＝その細胞について産生・受容を判断して書いた分子数
     （「なし」と判断した level 0 も含む）。<strong>未評価</strong>＝まだ手を付けていない分子数。<br>
     産生・受容の数字は level 1 以上、つまり実際に辺がある分子の数。
     <strong>辺の少なさは通信の少なさを意味しない</strong>——未評価が多ければ、単に手が回っていないだけ。
     薄く表示した行は評価が半分に届いていない細胞。</p>
  <div class="scroll-x">
  <table class="rel">
    <thead><tr><th>細胞</th><th>系</th><th class="num">評価済み</th><th class="num">産生</th>
      <th class="num">受容</th><th class="num">未評価</th></tr></thead>
    <tbody>{{cov_rows}}</tbody>
  </table>
  </div>
</section>""".replace("{cov_rows}", cov_rows)

    todos = [r for r in rows if r["status"] == "todo"]
    infs = [r for r in rows if r["status"] == "inferred"]
    figs = [r for r in rows if r["status"] == "figure_read"]

    def fig_table(rs):
        if not rs:
            return ('<p class="empty">なし。図から読み取った記述は現時点で0件。</p>')
        return ("<div class='scroll-x'><table class='rel todo'><thead><tr>"
                "<th>種別</th><th>対象</th><th>位置</th><th>読んだ図</th>"
                "<th>内容</th><th>出典</th></tr></thead><tbody>" + "".join(
                    f'<tr><td>{e(r["kind"])}</td>'
                    f'<td><a href="{ctx.url(r["kind"], r["id"])}">{e(r["id"])}</a></td>'
                    f'<td><code>{e(r["path"])}</code></td>'
                    f'<td><b>{e(r["figure"])}</b></td>'
                    f'<td>{e(r["label"])}</td>'
                    f'<td>{src_link(r["source"])}</td></tr>' for r in rs)
                + "</tbody></table></div>")
    return f"""
<h1>TODO</h1>
<p class="lede">未検証がどこに残っているかが常に見える状態を保つ。
   これが崩れた時点でこのサイトは信用できなくなる。</p>

<section class="card">
  <div class="counts">
    <span class="cnt verified">verified {counts.get("verified", 0)}</span>
    <span class="cnt figure">図から読取 {counts.get("figure_read", 0)}</span>
    <span class="cnt inferred">inferred {counts.get("inferred", 0)}</span>
    <span class="cnt todo">todo {counts.get("todo", 0)}</span>
  </div>
</section>

{cov_html}

<section class="card">
  <h2>図から読み取ったもの <span class="sub">{len(figs)}</span></h2>
  <p class="note"><strong>本文に明記がなく、図を目視で読み取った記述。</strong>
     出典はあるが、根拠は私の視覚的解釈であって引用できる一文ではない。
     verified と同じ扱いにすると精度の違いが消えるので分けている。
     どの図を読んだかを必ず記録させており（無いとビルドが止まる）、後から検算できる。
     <strong>プライマー配列と抗体は、この status でもサイトに出力しない</strong>
     （出力は verified のみ）。</p>
  {fig_table(figs)}
</section>

<section class="card">
  <h2>status: todo <span class="sub">{len(todos)}</span></h2>
  <p class="note">未着手・要調査。空欄のまま出力されている。</p>
  {table(todos)}
</section>

<section class="card">
  <h2>status: inferred <span class="sub">{len(infs)}</span></h2>
  <p class="note">一般知識から書かれたもの。出典を当てて verified に上げるか、消す。</p>
  {table(infs)}
</section>"""


# --------------------------------------------------------------------------
CSS = """
:root{--fg:#1b1f24;--muted:#6b7480;--line:#e2e6ea;--bg:#fbfcfd;--card:#fff;
--accent:#2b6cb0;--warn:#b45309;--alarm:#b91c1c;--ok:#166534}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.75 -apple-system,BlinkMacSystemFont,"Hiragino Sans","Noto Sans JP",sans-serif}
main{max-width:1040px;margin:0 auto;padding:24px 20px 64px}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
h1{font-size:26px;margin:8px 0 4px;letter-spacing:.01em}
h2{font-size:17px;margin:0 0 12px;padding-bottom:8px;border-bottom:1px solid var(--line)}
h3{font-size:14px;margin:20px 0 8px;color:var(--muted)}
.sub{font-weight:400;font-size:12px;color:var(--muted)}
.lede{color:var(--muted);margin:0 0 20px}
.note{font-size:13px;color:var(--muted)}
.empty{color:var(--muted);font-size:13px;font-style:italic}
code{background:#f1f3f5;padding:1px 5px;border-radius:3px;font-size:12px}
.topbar{display:flex;align-items:center;gap:18px;padding:12px 20px;background:var(--card);
border-bottom:1px solid var(--line);flex-wrap:wrap}
.brand{font-weight:700;color:var(--fg)}
.topbar nav{display:flex;gap:14px}
.topbar nav a{font-size:14px;color:var(--muted)}
.topbar nav a.active{color:var(--accent);font-weight:600}
.flag{margin-left:auto;font-size:11px;color:var(--warn);border:1px solid currentColor;
border-radius:3px;padding:1px 7px}
footer{border-top:1px solid var(--line);padding:20px;text-align:center;
font-size:12px;color:var(--muted)}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;
padding:18px 20px;margin:0 0 18px}
.scroll-x{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13px}
.rel th,.rel td{border-bottom:1px solid var(--line);padding:7px 10px;text-align:left;
vertical-align:top}
.rel th{color:var(--muted);font-weight:600;font-size:12px}
.rel .num{text-align:right;font-variant-numeric:tabular-nums}
/* 行列 */
.matrix{font-size:12px;width:auto}   /* table{width:100%} を打ち消す。伸ばすと行見出しが間延びする */
.matrix th{padding:6px 9px;color:var(--muted);font-weight:600;white-space:nowrap}
.matrix tbody th{text-align:right}
.matrix td{border:1px solid var(--line);width:64px;height:40px;text-align:center}
.matrix td.m0{background:repeating-linear-gradient(45deg,#fafbfc,#fafbfc 4px,#f2f4f6 4px,#f2f4f6 8px)}
.matrix td.mx{cursor:pointer}
.matrix td.mx:hover,.matrix td.mx:focus{outline:2px solid var(--accent)}
.matrix td.sel{outline:2px solid var(--accent)}
.h1{background:#e8f1fb}.h2{background:#d3e5f8}.h3{background:#bad6f4}
.h4{background:#9fc6ef}.h5{background:#82b4e9}
.matrix .n{font-weight:600}
.matrix-status{font-size:13px;min-height:1.5em}
.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:5px;
vertical-align:middle}
.dot.big{width:14px;height:14px}
/* 絞り込み */
.filters{display:flex;gap:14px;flex-wrap:wrap;align-items:flex-end}
.filters label{display:flex;flex-direction:column;font-size:12px;color:var(--muted);gap:4px}
.filters input,.filters select{font:inherit;font-size:13px;padding:5px 8px;
border:1px solid var(--line);border-radius:5px;background:#fff;min-width:150px}
.filters button{font:inherit;font-size:13px;padding:6px 14px;border:1px solid var(--line);
border-radius:5px;background:#fff;cursor:pointer}
/* チップ・細胞 */
.gene{font-size:11px;color:var(--muted);font-weight:400;margin-left:6px}
.chips{display:flex;flex-wrap:wrap;gap:5px;margin:6px 0}
.chip{font-size:11px;background:#eef2f6;border-radius:11px;padding:2px 9px;color:#41505f}
.chip.ghost{opacity:.45}
.cellchip{display:inline-flex;align-items:center;font-size:13px}
.outside{font-size:10px;color:var(--warn);border:1px solid currentColor;border-radius:3px;
padding:0 4px;margin-left:5px}
/* 閉じ方 */
.closure{display:inline-block;font-size:11px;border-radius:3px;padding:1px 8px;margin:0}
.closure.big{font-size:13px;padding:4px 12px;margin:6px 0 12px}
.c-closed{background:#e7f5ec;color:var(--ok)}
.c-out_only{background:#fdf0e3;color:var(--warn)}
.c-in_only{background:#eef2f6;color:#41505f}
.c-external{background:#f1f3f5;color:var(--muted)}
.alarm{color:var(--alarm);font-size:14px;border-left:3px solid currentColor;
padding:6px 12px;background:#fef2f2;margin:0 0 16px}
/* level / effects */
.lv{display:inline-block;width:18px;height:18px;line-height:18px;text-align:center;
border-radius:3px;font-size:11px;font-weight:700;color:#fff}
.lv3{background:#1d4ed8}.lv2{background:#60a5fa}.lv1{background:#bfdbfe;color:#1e3a5f}
.lv0{background:#e5e7eb;color:var(--muted)}
/* 効果の色分け: 負=青 / ±0=灰 / 正=赤。濃さが絶対値に対応する。
   クラス名に '+' は使わない（CSS セレクタでエスケープが要るため）。 */
.effv{display:inline-block;min-width:34px;padding:1px 6px;border-radius:3px;
font-weight:700;font-variant-numeric:tabular-nums;text-align:center}
.effv.vm2{background:#1d4ed8;color:#fff}
.effv.vm1{background:#dbeafe;color:#1e40af}
.effv.v0{background:transparent;color:#b6bec7;font-weight:400}
.effv.vp1{background:#fee2e2;color:#b91c1c}
.effv.vp2{background:#b91c1c;color:#fff}
.eff-table td{padding:3px 8px;border:0}
.effcell{display:flex;align-items:center;gap:8px}
.bar{display:inline-block;height:8px;border-radius:2px}
.bar.vm2{width:48px;background:#1d4ed8}.bar.vm1{width:24px;background:#93c5fd}
.bar.v0{width:6px;background:#d1d5db}.bar.vp1{width:24px;background:#fca5a5}
.bar.vp2{width:48px;background:#b91c1c}
/* 分子一覧テーブル
   見出しを固定するには、縦スクロールがこの箱の中で起きている必要がある。
   ページ側でスクロールすると sticky は overflow 親の中で効かず、見出しが流れてしまう。 */
.mol-scroll{max-height:70vh;overflow:auto;border:1px solid var(--line);border-radius:6px}
.mol-table{font-size:13px;width:100%}
.mol-table th,.mol-table td{border-bottom:1px solid var(--line);padding:6px 8px;
text-align:left;white-space:nowrap}
.mol-table td .na{color:#cfd6dd}
.mol-table thead th{color:var(--muted);font-size:12px;font-weight:600;
border-bottom:2px solid var(--line);position:sticky;top:var(--legend-h,0px);z-index:2;
background:var(--card);box-shadow:0 1px 0 var(--line)}
/* 分子名の列は横スクロール時に左端へ固定する */
.mol-table td.sym{position:sticky;left:0;z-index:1;background:var(--card)}
.mol-table tbody tr:hover td.sym{background:#f7fafd}
.mol-table thead th:first-child{left:0;z-index:3}
.mol-table th.sortable{cursor:pointer;user-select:none}
.mol-table th.sortable:hover{color:var(--accent)}
.mol-table th.num,.mol-table td.num{text-align:right}
.mol-table .arrow{display:inline-block;width:12px;color:var(--accent)}
.mol-table th.asc .arrow::after{content:"▲";font-size:9px}
.mol-table th.desc .arrow::after{content:"▼";font-size:9px}
.mol-table tbody tr:hover{background:#f7fafd}
.mol-table .sym a{font-weight:600}
.mol-table .gene-col{color:var(--muted);font-size:12px}
.mol-table .dots{margin-left:6px;white-space:nowrap}
.mol-table .dots .dot{margin-right:1px}
.origin{display:block;font-size:10px;font-weight:400;color:var(--muted);
margin-top:1px;white-space:nowrap}
/* 色と細胞の対応表。スクロール領域の中に置き、見出しと一緒に上端へ固定する。
   外に置くと、画面の低いノートPCではページを送ったときに凡例だけ視界から消える。 */
.legend{position:sticky;top:0;z-index:4;border-bottom:1px solid var(--line);
padding:8px 12px;background:#fcfdfe}
.legend-row{display:flex;flex-wrap:wrap;align-items:center;gap:4px 10px;padding:3px 0}
.legend-head{font-size:11px;font-weight:700;color:var(--muted);min-width:34px}
.legend-head.out{color:var(--warn)}
.legend-chip{display:inline-flex;align-items:baseline;gap:3px;font-size:12px;color:var(--fg);
padding:1px 6px;border-radius:11px}
.legend-chip:hover{background:#eef2f6;text-decoration:none}
.legend-origin{font-size:10px;color:var(--muted)}
.legend-note{font-size:10px;color:var(--muted);flex-basis:100%}
.disclaimer{border:1px solid #d9b28c;background:#fffaf3;border-radius:7px;
padding:12px 16px;margin:0 0 18px;font-size:13px;line-height:1.9}
.unwritten{font-size:13px;color:var(--warn);background:#fffaf3;border-left:3px solid currentColor;
padding:8px 12px;border-radius:0 4px 4px 0}
.assessed-zero{font-size:12px;color:var(--muted);margin-top:10px;padding-top:8px;
border-top:1px dashed var(--line)}
.unassessed-n{color:var(--warn)}
/* 類洞の断面図 */
.sn-fig{margin:0 0 16px}
.sn-fig figcaption{font-size:12px;color:var(--muted);text-align:center;margin-top:8px;
line-height:1.7}
.sn-svg{width:100%;max-width:760px;height:auto;display:block;margin:0 auto;color:#7c8792}
.sn-svg .sn-bg{fill:#fdfdfe}
.sn-svg .sn-hep{fill:#fdf3e9}
.sn-svg .sn-elsewhere{fill:#f2f4f6}
.sn-svg .sn-disse{fill:#f6f2fb}
.sn-svg .sn-lumen{fill:#eaf4fb}
.sn-svg .sn-hepcell{fill:#f7e3cd;stroke:#e0b98d;stroke-width:1.2}
.sn-svg .sn-endo{stroke:#4CC9F0;stroke-width:3;stroke-dasharray:14 7;stroke-linecap:round}
.sn-svg .sn-zone2{font:10px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#a8b2bd}
.sn-svg .sn-zone{font:600 11px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#8b95a1;
letter-spacing:.06em}
.sn-svg .sn-cell{stroke-width:1.8;fill:#fff;stroke:#9aa5b1}
.sn-svg .sn-cell.r-prod{fill:#dcefff;stroke:#3f7fbf}
.sn-svg .sn-cell.r-recv{fill:#ffe6e6;stroke:#c96a6a}
.sn-svg .sn-cell.r-both{fill:#e6f6ea;stroke:#3f9c63}
.sn-svg .sn-cell.sn-outside{stroke-dasharray:5 4}
.sn-svg .sn-cell.is-focus{stroke-width:3.4;stroke:#1d4ed8}
.sn-svg .sn-name{font:600 12px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#1b1f24}
.sn-svg .sn-sub{font:10px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#7c8792}
.sn-svg a:hover .sn-name{fill:var(--accent);text-decoration:underline}
.sn-svg .sn-edge{fill:none;stroke-width:1.5;opacity:.8}
.sn-svg .sn-edge.out{stroke:#3f7fbf}
.sn-svg .sn-edge.in{stroke:#c96a6a}
.sn-svg .sn-edge.lv3{stroke-width:3.2;opacity:1}
.sn-svg .sn-edge.lv2{stroke-width:2.2;opacity:.9}
.sn-svg .sn-hub{fill:#fff;stroke:#41505f;stroke-width:2}
.sn-svg .sn-hubname{font:700 13px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#1b1f24}
/* 分子の形（会合状態） */
.oli{height:26px;width:auto;vertical-align:middle;color:#5b8fc9}
.oli-line{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px;
font-size:13px}
.oli-u0{fill:#bcd8f2;stroke:#4a7fb5;stroke-width:1.4}
.oli-u1{fill:#f6dcc0;stroke:#c08a44;stroke-width:1.4}
.pw-svg .oli-u0{fill:#bcd8f2;stroke:#4a7fb5;stroke-width:1.2}
.pw-svg .oli-u1{fill:#f6dcc0;stroke:#c08a44;stroke-width:1.2}
/* 経路図: 細胞の解剖図として描く */
.pw-fig{margin:0}
.pw-fig figcaption{font-size:12px;color:var(--muted);text-align:center;margin-top:10px;
line-height:1.7}
.pw-svg{width:100%;max-width:760px;height:auto;display:block;margin:0 auto;color:#6b7480}
.pw-svg .pw-outside{fill:#f4f8fc}
.pw-svg .pw-cytosol{fill:#fdfbf6}
.pw-svg .pw-mem{stroke:#c9b28a;stroke-width:1.6}
.pw-svg .pw-lipid{fill:#e0cba6}
.pw-svg .pw-nucleus{fill:#f3eefb;stroke:#c3b2dd;stroke-width:1.4}
.pw-svg .pw-zone{font:600 11px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#9aa5b1;
letter-spacing:.08em}
.pw-svg .pw-box{stroke-width:1.5}
.pw-svg .pw-box.k-ligand{fill:#e8f1fb;stroke:#5b8fc9}
.pw-svg .pw-box.k-receptor{fill:#e7f5ec;stroke:#4aa373}
.pw-svg .pw-box.k-transducer{fill:#f6f1fd;stroke:#8e6fc0}
.pw-svg .pw-box.k-tf{fill:#fdf0e3;stroke:#c98a3d}
.pw-svg .pw-box.k-output{fill:#fdeaea;stroke:#c96a6a}
.pw-svg .pw-kind{font:600 10px sans-serif;fill:#6b7480;letter-spacing:.04em}
.pw-svg .pw-label{font:600 13px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#1b1f24}
.pw-svg .pw-mol{font:11px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#5c6672}
.pw-svg .pw-arrow{stroke:currentColor;stroke-width:1.5;fill:none}
/* 受容体の膜での形 */
.pw-svg .rc-chain{fill:#4aa373;stroke:#2f7551;stroke-width:1.2}
.pw-svg .rc-chain2{fill:#bfe3cd;stroke:#4aa373;stroke-width:1.2}
.pw-svg .rc-loop{fill:none;stroke:#2f7551;stroke-width:2.4}
.pw-svg .rc-kin{fill:#e8dcf8;stroke:#8e6fc0;stroke-width:1.2}
.pw-svg .rc-jak{fill:#e8dcf8;stroke:#8e6fc0;stroke-width:1.2}
.pw-svg .rc-g{fill:#fdf0e3;stroke:#c98a3d;stroke-width:1.2}
.pw-svg .rc-ig{fill:#dff0e7;stroke:#4aa373;stroke-width:1.1}
.pw-svg .rc-mini{font:600 8px sans-serif;fill:#41505f}
/* 産生経路: 受容側と逆向きなので区画の上下が入れ替わる */
.pd-svg .pd-cytosol{fill:#fdfbf6}
.pd-svg .pd-outside{fill:#f4f8fc}
.pd-svg .pd-nucleus{fill:#f3eefb;stroke:#c3b2dd;stroke-width:1.4}
.pd-svg .pd-zone{font:600 11px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#9aa5b1;
letter-spacing:.08em}
.pd-tag{display:inline-block;font-size:10px;border-radius:3px;padding:1px 7px;
white-space:nowrap;border:1px solid}
.p-trigger{background:#fdf0e3;border-color:#c98a3d;color:#94602a}
.p-transcription{background:#f3eefb;border-color:#8e6fc0;color:#654a94}
.p-proform{background:#eef2f6;border-color:#8b95a1;color:#41505f}
.p-processing{background:#fdeaea;border-color:#c96a6a;color:#a04747}
.p-release{background:#e8f1fb;border-color:#5b8fc9;color:#2c5f96}
.p-extracellular{background:#e7f5ec;border-color:#4aa373;color:#2f7551}
.arch-table td{vertical-align:top}
.arch-table td:nth-child(2),.arch-table td:nth-child(3){min-width:230px;font-size:12px;
line-height:1.8}
.arch-table tr:target td{background:#fffaf3}
.pw-svg .pw-edge{font:11px "Hiragino Sans","Noto Sans JP",sans-serif;fill:#5c6672;
paint-order:stroke;stroke:#fff;stroke-width:4px}
.pw-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:14px}
.pw-card{border:1px solid var(--line);border-radius:7px;padding:12px 14px;background:#fff}
.pw-card h3{margin:0 0 4px;font-size:14px;color:var(--fg)}
.pw-tag{display:inline-block;font-size:10px;border-radius:3px;padding:1px 7px;
white-space:nowrap;border:1px solid}
.k-ligand{background:#e8f1fb;border-color:#5b8fc9;color:#2c5f96}
.k-receptor{background:#e7f5ec;border-color:#4aa373;color:#2f7551}
.k-transducer{background:#f3eefb;border-color:#8e6fc0;color:#654a94}
.k-tf{background:#fdf0e3;border-color:#c98a3d;color:#94602a}
.k-output{background:#fdeaea;border-color:#c96a6a;color:#a04747}
.pw-steps td:first-child{width:70px;vertical-align:top}
.pw-steps p{margin:3px 0}
.outcome{margin:0 0 8px;font-size:14px;line-height:1.9}
/* 分子間の関係 */
.rel-arrow{font-size:17px;color:var(--accent);text-align:center;width:34px}
.rel-tag{font-size:11px;background:#eef2f6;border-radius:11px;padding:1px 9px;
white-space:nowrap;color:#41505f}
tr.cov-thin td{background:#fffaf3}
tr.cov-thin td:first-child{box-shadow:inset 3px 0 0 #d9a566}
/* status 表示（CLAUDE.md §0.2） */
.block{margin:0 0 10px;padding:2px 0}
.block.is-inferred{color:var(--muted);border:1px dotted #c3ccd5;border-radius:5px;
padding:8px 12px;background:#fcfcfd}
.block.is-todo{border:1px dashed #d9b28c;border-radius:5px;padding:8px 12px;background:#fffdf8}
tr.is-inferred td{color:var(--muted)}
.badge{font-size:10px;border-radius:3px;padding:1px 6px;margin-left:6px;white-space:nowrap}
.badge-inferred{background:#eef2f6;color:var(--muted);border:1px dotted #b9c4cf}
.badge-todo{background:#fff4e5;color:var(--warn);border:1px dashed currentColor}
.badge-figure{background:#eaf2fb;color:#2c5f96;border:1px solid #9dbfe0}
.block.is-figure{border:1px solid #cfe0f2;border-radius:5px;padding:8px 12px;background:#fbfdff}
.src{font-size:11px;background:#eef4fb;border-radius:3px;padding:1px 6px;white-space:nowrap}
/* ファミリー */
.tabs{display:flex;gap:6px;margin-bottom:16px;flex-wrap:wrap}
.tab{font:inherit;font-size:13px;padding:6px 14px;border:1px solid var(--line);
background:#fff;border-radius:5px;cursor:pointer}
.tab.active{background:var(--accent);color:#fff;border-color:var(--accent)}
.fam-node{border-left:3px solid var(--line);padding:8px 0 8px 14px;margin:0 0 6px}
.fam-node.d1{margin-left:26px}.fam-node.d2{margin-left:52px}
.fam-name{font-weight:600;font-size:14px}
.fam-node p{margin:2px 0 6px}
.mismatch td:nth-child(4){min-width:280px}
.own{border-color:#c7a87a;background:#fffdf7}
.counts{display:flex;gap:10px;flex-wrap:wrap}
.cnt{font-size:13px;border-radius:4px;padding:4px 12px}
.cnt.verified{background:#e7f5ec;color:var(--ok)}
.cnt.figure{background:#eaf2fb;color:#2c5f96}
.cnt.inferred{background:#eef2f6;color:var(--muted)}
.cnt.todo{background:#fff4e5;color:var(--warn)}
.refs{font-size:13px;padding-left:20px}
@media(max-width:640px){main{padding:16px 12px 48px}}
"""

JS = """
(function(){
  var payload = document.getElementById('atlas-data');
  var DATA = payload ? JSON.parse(payload.textContent) : null;

  /* --- ファミリーの軸切り替え --- */
  var tabs = document.querySelectorAll('.tab');
  if (tabs.length) {
    function show(axis){
      tabs.forEach(function(t){ t.classList.toggle('active', t.dataset.axis === axis); });
      document.querySelectorAll('.axis-panel').forEach(function(p){
        p.hidden = (p.dataset.axis !== axis);
      });
      history.replaceState(null,'','#'+axis);
    }
    tabs.forEach(function(t){ t.addEventListener('click', function(){ show(t.dataset.axis); }); });
    var initial = location.hash.replace('#','') || tabs[0].dataset.axis;
    show(document.querySelector('.tab[data-axis="'+initial+'"]') ? initial : tabs[0].dataset.axis);
  }

  /* --- 分子間の関係の絞り込み --- */
  var relSel = document.getElementById('rel-filter');
  if (relSel) {
    relSel.addEventListener('input', function(){
      var want = relSel.value, shown = 0;
      document.querySelectorAll('#rel-body tr.rel-row').forEach(function(tr){
        var ok = !want || tr.dataset.type === want;
        tr.hidden = !ok;
        if (ok) shown++;
      });
      document.getElementById('rel-count').textContent = shown;
    });
  }

  /* --- トップページの絞り込み --- */
  var grid = document.getElementById('mol-grid');
  if (!grid || !DATA) return;
  var state = {q:'', closure:'', status:'', pair:null};
  ['structure','receptor','naming'].forEach(function(ax){ state['fam-'+ax] = ''; });

  function apply(){
    var shown = 0;
    document.querySelectorAll('.mol').forEach(function(el){
      var d = DATA.molecules[el.dataset.id];
      var ok = true;
      if (state.q) {
        var hay = (d.search || '').toLowerCase();
        ok = ok && hay.indexOf(state.q.toLowerCase()) !== -1;
      }
      ['structure','receptor','naming'].forEach(function(ax){
        var want = state['fam-'+ax];
        if (want) ok = ok && (d.families || {})[ax] === want;
      });
      if (state.closure) ok = ok && d.closure === state.closure;
      if (state.status === 'verified')     ok = ok && d.counts.verified > 0;
      if (state.status === 'has-inferred') ok = ok && d.counts.inferred > 0;
      if (state.status === 'has-todo')     ok = ok && d.counts.todo > 0;
      if (state.pair) {
        var cell = ((DATA.matrix[state.pair[0]] || {})[state.pair[1]]) || [];
        ok = ok && cell.indexOf(el.dataset.id) !== -1;
      }
      el.hidden = !ok;
      if (ok) shown++;
    });
    var c = document.getElementById('count');
    if (c) c.textContent = shown + (shown === DATA.count ? '' : ' / ' + DATA.count);
  }

  document.querySelectorAll('[data-filter]').forEach(function(el){
    el.addEventListener('input', function(){ state[el.dataset.filter] = el.value; apply(); });
  });

  var status = document.getElementById('matrix-status');
  document.querySelectorAll('.matrix td.mx').forEach(function(td){
    function pick(){
      var same = state.pair && state.pair[0] === td.dataset.producer
                            && state.pair[1] === td.dataset.receiver;
      document.querySelectorAll('.matrix td').forEach(function(x){ x.classList.remove('sel'); });
      if (same) { state.pair = null; status.textContent = ''; }
      else {
        state.pair = [td.dataset.producer, td.dataset.receiver];
        td.classList.add('sel');
        status.textContent = DATA.cellNames[state.pair[0]] + ' → ' +
          DATA.cellNames[state.pair[1]] + ' の ' + td.dataset.mols.split(',').length +
          ' 分子に絞り込み中（もう一度クリックで解除）';
      }
      apply();
      grid.scrollIntoView({behavior:'smooth', block:'nearest'});
    }
    td.addEventListener('click', pick);
    td.addEventListener('keydown', function(ev){
      if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); pick(); }
    });
  });

  /* --- 凡例の高さ分だけ見出しの固定位置を下げる --- */
  var legend = document.querySelector('.mol-scroll .legend');
  var scrollBox = document.querySelector('.mol-scroll');
  if (legend && scrollBox) {
    var setLegendH = function(){
      scrollBox.style.setProperty('--legend-h', legend.offsetHeight + 'px');
    };
    setLegendH();
    window.addEventListener('resize', setLegendH);
    if (window.ResizeObserver) new ResizeObserver(setLegendH).observe(legend);
  }

  /* --- 列見出しクリックで並び替え --- */
  var table = document.getElementById('mol-table');
  if (table) {
    var tbody = document.getElementById('mol-grid');
    var origOrder = [].slice.call(tbody.querySelectorAll('tr'));
    table.querySelectorAll('th.sortable').forEach(function(th){
      function sortBy(){
        var col = +th.dataset.col, numeric = th.dataset.type === 'num';
        var wasAsc = th.classList.contains('asc');
        var wasDesc = th.classList.contains('desc');
        table.querySelectorAll('th').forEach(function(x){
          x.classList.remove('asc','desc');
        });
        if (wasDesc) {                       /* 3回目のクリックで元の順に戻す */
          origOrder.forEach(function(tr){ tbody.appendChild(tr); });
          return;
        }
        var asc = !wasAsc;
        th.classList.add(asc ? 'asc' : 'desc');
        var rows = [].slice.call(tbody.querySelectorAll('tr'));
        rows.sort(function(a, b){
          var ca = a.children[col], cb = b.children[col];
          var va, vb;
          if (numeric) {
            va = parseFloat(ca.dataset.v); vb = parseFloat(cb.dataset.v);
            if (isNaN(va)) va = -Infinity;
            if (isNaN(vb)) vb = -Infinity;
          } else {
            va = ca.textContent.trim(); vb = cb.textContent.trim();
            return (asc ? 1 : -1) * va.localeCompare(vb, 'ja');
          }
          if (va === vb) {                   /* 同値は分子名で安定させる */
            return a.children[0].textContent.localeCompare(
                   b.children[0].textContent, 'ja');
          }
          return (asc ? 1 : -1) * (va - vb);
        });
        rows.forEach(function(tr){ tbody.appendChild(tr); });
      }
      th.addEventListener('click', sortBy);
      th.addEventListener('keydown', function(ev){
        if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); sortBy(); }
      });
    });
  }

  var reset = document.getElementById('reset');
  if (reset) reset.addEventListener('click', function(){
    document.querySelectorAll('[data-filter]').forEach(function(el){ el.value = ''; });
    Object.keys(state).forEach(function(k){ state[k] = (k === 'pair') ? null : ''; });
    document.querySelectorAll('.matrix td').forEach(function(x){ x.classList.remove('sel'); });
    status.textContent = '';
    apply();
  });

  apply();
})();
"""


# --------------------------------------------------------------------------
def build(skip_validate: bool) -> int:
    if not skip_validate:
        issues, _ = run_validate(strict=False, check_refs=False)
        errs = [i for i in issues if i.level == "error"]
        if errs:
            print(f"検証で ERROR {len(errs)} 件。ビルドしない。")
            for i in errs:
                print(i.fmt())
            print("\npython3 scripts/validate.py で全件を確認する。")
            return 1

    issues: list = []
    ds = load_dataset(issues)
    cfg = ds.config
    include_own = bool(cfg.get("include_own_data"))
    # GitHub Pages のプロジェクトページは /<repo>/ 配下に出る。
    # ワークフローが configure-pages の base_path を BASE_URL で渡してくる。
    base = os.environ.get("BASE_URL") or (cfg.get("site") or {}).get("base_url") or ""
    ctx = Ctx(base)
    inside = in_system_set(ds)

    mols = {mid: public_molecule(m, include_own) for mid, (m, _) in ds.molecules.items()}
    closures = {mid: derive_closure(m, inside) for mid, m in mols.items()}

    if SITE.exists():
        shutil.rmtree(SITE)
    SITE.mkdir(parents=True)
    write(SITE / "assets" / "style.css", CSS)
    write(SITE / "assets" / "app.js", JS)

    # --- クライアント用データ ---
    counts_all = {"verified": 0, "figure_read": 0, "inferred": 0, "todo": 0}
    payload_mols = {}
    for mid, m in mols.items():
        blocks: list = []
        # 件数は work_view（未検証のプライマー・抗体も数える）で数える。
        # 「TODO を含む」で絞ったときに、出力から落ちた作業が漏れないようにするため。
        walk_status_blocks(work_view(ds.molecules[mid][0], include_own), "", blocks, set())
        c = {"verified": 0, "figure_read": 0, "inferred": 0, "todo": 0}
        for _, w in blocks:
            st = w["_node"].get("status")
            if st in c:
                c[st] += 1
                counts_all[st] += 1
        search = " ".join(filter(None, [
            mid, m.get("symbol"), m.get("gene"), *(m.get("aliases") or []),
            *[ds.families[f][0].get("label_ja", "")
              for f in (m.get("families") or {}).values() if f in ds.families],
        ]))
        payload_mols[mid] = {"families": m.get("families") or {},
                             "closure": closures[mid]["key"],
                             "counts": c, "search": search}
    for store in (ds.cells, ds.families):
        for _, (obj, _f) in store.items():
            blocks = []
            walk_status_blocks(obj, "", blocks, set())
            for _, w in blocks:
                st = w["_node"].get("status")
                if st in counts_all:
                    counts_all[st] += 1

    cell_order = (cfg.get("matrix") or {}).get("cells") or []
    payload = {
        "count": len(mols),
        "molecules": payload_mols,
        "matrix": build_matrix(ds, mols, cell_order,
                               (cfg.get("matrix") or {}).get("min_level", 1)),
        "cellNames": {cid: f"{cell_generic(ds, cid)}（{cell_detail(ds, cid)}）"
                      for cid in ds.cells},
        "include_own_data": include_own,
    }
    write(SITE / "data.json", json.dumps(payload, ensure_ascii=False, indent=1))
    inline = ('<script type="application/json" id="atlas-data">'
              + json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")
              + "</script>")

    # --- ページ ---
    write(SITE / "index.html",
          page(ctx, "分子一覧",
               render_index(ctx, ds, mols, closures, cfg) + inline, cfg, ""))
    for mid, m in mols.items():
        write(SITE / "molecule" / mid / "index.html",
              page(ctx, m.get("symbol") or mid,
                   render_molecule(ctx, ds, mid, m, closures[mid], cfg), cfg))
    for cid, (c, _) in ds.cells.items():
        write(SITE / "cell" / cid / "index.html",
              page(ctx, f"{cell_generic(ds, cid)}（{cell_detail(ds, cid)}）",
                   render_cell(ctx, ds, cid, c, mols, inside), cfg))
    write(SITE / "pathway" / "index.html",
          page(ctx, "経路", render_pathway_index(ctx, ds, mols), cfg, "pathway"))
    for pid, (pw, _) in ds.pathways.items():
        write(SITE / "pathway" / pid / "index.html",
              page(ctx, pw.get("label_ja") or pid,
                   render_pathway(ctx, ds, pid, pw, mols), cfg, "pathway"))
    write(SITE / "production" / "index.html",
          page(ctx, "産生経路", render_production_index(ctx, ds, mols), cfg, "production"))
    for pid, (pr, _) in ds.production.items():
        write(SITE / "production" / pid / "index.html",
              page(ctx, pr.get("label_ja") or pid,
                   render_production(ctx, ds, pid, pr, mols), cfg, "production"))
    write(SITE / "architecture" / "index.html",
          page(ctx, "受容体の形", render_architectures(ctx, ds), cfg, "architecture"))
    write(SITE / "relations" / "index.html",
          page(ctx, "分子どうしの関係", render_relations(ctx, ds, mols), cfg, "relations"))
    write(SITE / "families" / "index.html",
          page(ctx, "ファミリー", render_families(ctx, ds, mols), cfg, "families"))
    for fid, (f, _) in ds.families.items():
        write(SITE / "family" / fid / "index.html",
              page(ctx, f.get("label_ja") or fid,
                   render_family(ctx, ds, fid, f, mols), cfg))
    todo_rows = collect_todo(ds, include_own)
    write(SITE / "todo" / "index.html",
          page(ctx, "TODO",
               render_todo(ctx, ds, todo_rows, counts_all, cell_coverage(ds, mols)),
               cfg, "todo"))

    n = sum(1 for _ in SITE.rglob("*.html"))
    print(f"site/ を生成: HTML {n} ページ"
          f"（分子 {len(mols)} / 細胞 {len(ds.cells)} / ファミリー {len(ds.families)}）")
    print(f"  verified {counts_all['verified']} / inferred {counts_all['inferred']} "
          f"/ todo {counts_all['todo']}")
    if not include_own:
        print("  include_own_data: false のため own_data は出力していない")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="data/ → site/")
    ap.add_argument("--skip-validate", action="store_true")
    ap.add_argument("--serve", action="store_true", help="生成後にローカルサーバを立てる")
    ap.add_argument("--port", type=int, default=8000)
    args = ap.parse_args()

    rc = build(args.skip_validate)
    if rc or not args.serve:
        return rc

    import functools
    import http.server
    import socketserver
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(SITE))
    with socketserver.TCPServer(("127.0.0.1", args.port), handler) as httpd:
        print(f"http://127.0.0.1:{args.port}/  (Ctrl-C で終了)")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
