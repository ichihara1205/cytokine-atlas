#!/usr/bin/env python3
"""data/ の整合性を検査する。CLAUDE.md §0.3 のビルド前チェック。

ERROR が1件でもあればビルドを止める（exit 1）。WARN は止めない。

    python3 scripts/validate.py
    python3 scripts/validate.py --strict        # WARN もエラー扱い
    python3 scripts/validate.py --check-refs    # PMID/DOI が実在するか PubMed/Crossref に問い合わせる（要ネットワーク）
    python3 scripts/validate.py --json          # 機械可読出力

ERROR にする条件（CLAUDE.md §0.3）:
  - status: verified / figure_read なのに source が空
  - status: figure_read なのに figure（どの図か）が空
  - 存在しない cell_id / pathway_id / theme_id / family_id を参照している
  - gene が HGNC 承認シンボルと一致しない
  - プライマー配列に ATGC 以外の文字が入っている
加えて、放置すると静かに壊れるものを同じ強度で止める:
  - id の重複 / ファイル名と id の不一致 / status 語彙外
  - system_closure の手書き（build.py が自動判定するので手で書くと必ず食い違う）
  - level・effects の範囲外 / 分化グラフの循環
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover
    sys.exit("PyYAML が必要です:  python3 -m pip install pyyaml")

ROOT = Path(__file__).resolve().parent.parent

# status の語彙（CLAUDE.md §0.2）。
# figure_read は「文献の図から読み取った」。本文に明記がなく目視解釈に依存するため、
# verified とは別扱いにする。出典に加えて、どの図かを figure に必ず書かせる。
STATUSES = {"verified", "figure_read", "inferred", "todo"}
# 出典が必須の status
NEEDS_SOURCE = {"verified", "figure_read"}
AXES = {"structure", "receptor", "naming"}
CELL_STAGES = {"pluripotent", "progenitor", "mature"}
CELL_ORIGINS = {"primary", "line", "ipsc_derived"}   # null = 培養していない参照ノード
EFFECT_KEYS = {"fibrosis", "inflammation", "lipid", "differentiation"}
SECRETION = {"soluble", "membrane_bound", "both", "latent"}
# 会合状態。図で形を描き分けるために持つ（単量体か二量体かで受容体の掴み方が変わる）
OLIGOMER = {"monomer", "homodimer", "homotrimer", "heterodimer", "multimer"}
# 受容体の膜での形は data/architectures/ に置く。
# 「その形でなければならない理由」は知識なので、status と出典を持てる data 側に住む。
TM_COUNTS = {1, 7, "multi"}
# 受容側の経路のステップ種別（外 → 中）
STEP_KINDS = {"ligand", "receptor", "transducer", "tf", "output"}
# 産生側の経路のステップ種別（中 → 外）。受容側とは向きも語彙も違うので分けて持つ
PROD_STEP_KINDS = {"trigger", "transcription", "proform", "processing",
                   "release", "extracellular"}
# 放出の手段
RELEASE_ROUTES = {"classical", "shedding", "gasdermin_pore", "cell_death", "unconventional"}
# 分子間辺の語彙。閉じた集合にしておかないと表記ゆれで同じ関係が別物になる
RELATION_TYPES = {
    "induces":          "下流で転写を誘導する",
    "antagonizes":      "受容体を奪うがシグナルを出さない",
    "sequesters":       "結合して遊離型を減らす",
    "heterodimer_with": "共有結合二量体としてのみ働く",
    "synergizes_with":  "単独では弱く、組むと強い",
    "opposes":          "作用が逆向き",
}
RANGE_MODES = {"autocrine", "juxtacrine", "paracrine", "endocrine"}

# own_data は本人の実測。外部 source は存在し得ないので verified でも source を要求しない。
# 代わりに date と実験の記述を要求する（下の check_own_data）。
NO_SOURCE_REQUIRED = {"own_data"}

PMID_RE = re.compile(r"^PMID:\s*(\d{1,9})$")
DOI_RE = re.compile(r"^(DOI:\s*)?10\.\d{4,9}/\S+$", re.I)


# --------------------------------------------------------------------------
# 収集
# --------------------------------------------------------------------------
@dataclass
class Issue:
    level: str          # error | warn
    file: str
    path: str           # YAML 内の位置
    msg: str
    hint: str | None = None

    def fmt(self) -> str:
        tag = "ERROR" if self.level == "error" else "WARN "
        loc = f"{self.file}:{self.path}" if self.path else self.file
        out = f"  [{tag}] {loc}\n          {self.msg}"
        if self.hint:
            out += f"\n          → {self.hint}"
        return out


@dataclass
class Dataset:
    config: dict = field(default_factory=dict)
    molecules: dict = field(default_factory=dict)   # id -> (data, relpath)
    cells: dict = field(default_factory=dict)
    families: dict = field(default_factory=dict)
    pathways: dict = field(default_factory=dict)
    themes: dict = field(default_factory=dict)
    architectures: dict = field(default_factory=dict)
    production: dict = field(default_factory=dict)


def _load_yaml(p: Path, issues: list[Issue]):
    try:
        return yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        issues.append(Issue("error", _rel(p), "", f"YAML として読めない: {e}"))
        return None


def _rel(p: Path) -> str:
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def _entries(doc, p: Path, issues: list[Issue]):
    """1ファイル=1件（マッピング）でも 1ファイル=複数件（リスト）でも受ける。"""
    if doc is None:
        return []
    if isinstance(doc, dict):
        return [doc]
    if isinstance(doc, list):
        out = []
        for i, item in enumerate(doc):
            if isinstance(item, dict):
                out.append(item)
            else:
                issues.append(Issue("error", _rel(p), f"[{i}]", "リストの要素がマッピングでない"))
        return out
    issues.append(Issue("error", _rel(p), "", "トップレベルはマッピングかリストであること"))
    return []


def load_dataset(issues: list[Issue]) -> Dataset:
    ds = Dataset()

    cfg = ROOT / "config.yaml"
    if cfg.exists():
        ds.config = _load_yaml(cfg, issues) or {}
    else:
        issues.append(Issue("error", "config.yaml", "", "config.yaml が無い"))

    for kind, attr in (
        ("molecules", "molecules"),
        ("cells", "cells"),
        ("families", "families"),
        ("pathways", "pathways"),
        ("themes", "themes"),
        ("architectures", "architectures"),
        ("production", "production"),
    ):
        d = ROOT / "data" / kind
        if not d.is_dir():
            continue
        store = getattr(ds, attr)
        for p in sorted(d.rglob("*.yaml")):
            if p.name.startswith("_"):      # _TEMPLATE.yaml などは読み込まない
                continue
            doc = _load_yaml(p, issues)
            items = _entries(doc, p, issues)
            for item in items:
                oid = item.get("id")
                if not oid:
                    issues.append(Issue("error", _rel(p), "id", "id が無い"))
                    continue
                if oid in store:
                    issues.append(Issue(
                        "error", _rel(p), f"id={oid}",
                        f"id が重複している（既出: {store[oid][1]}）"))
                    continue
                # 1ファイル1件のときはファイル名と id を一致させる
                if len(items) == 1 and p.stem != str(oid):
                    issues.append(Issue(
                        "warn", _rel(p), f"id={oid}",
                        f"ファイル名 '{p.stem}' と id '{oid}' が違う",
                        f"{oid}.yaml にリネームすると探しやすい"))
                store[str(oid)] = (item, _rel(p))
    return ds


def load_hgnc(ds: Dataset, issues: list[Issue]) -> dict[str, tuple[str, str]]:
    """symbol -> (kind, current_symbol)"""
    rel = (ds.config.get("validate") or {}).get("hgnc_file") or "data/hgnc/symbols.tsv"
    p = ROOT / rel
    if not p.exists():
        issues.append(Issue(
            "error", rel, "",
            "HGNC 照合表が無いので gene のチェックができない",
            "python3 scripts/fetch_hgnc.py で生成する"))
        return {}
    table: dict[str, tuple[str, str]] = {}
    for line in p.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) == 3:
            table[parts[0]] = (parts[1], parts[2])
    return table


# --------------------------------------------------------------------------
# 汎用チェック: status / source
# --------------------------------------------------------------------------
def walk_status_blocks(node, path: str, out: list[tuple[str, dict]], under: set[str]):
    """status キーを持つ辞書を再帰的に集める。under は祖先のキー名の集合。"""
    if isinstance(node, dict):
        if "status" in node:
            out.append((path, {"_node": node, "_under": set(under)}))
        for k, v in node.items():
            walk_status_blocks(v, f"{path}.{k}" if path else str(k), out, under | {str(k)})
    elif isinstance(node, list):
        for i, v in enumerate(node):
            walk_status_blocks(v, f"{path}[{i}]", out, under)


def _empty(v) -> bool:
    return v is None or (isinstance(v, str) and not v.strip()) or v == []


def check_status_and_sources(obj: dict, relfile: str, issues: list[Issue], strict_inferred: bool):
    blocks: list[tuple[str, dict]] = []
    walk_status_blocks(obj, "", blocks, set())
    for path, wrapped in blocks:
        node, under = wrapped["_node"], wrapped["_under"]
        st = node.get("status")
        if st not in STATUSES:
            issues.append(Issue(
                "error", relfile, path or "status",
                f"status: {st!r} は不正",
                f"{' / '.join(sorted(STATUSES))} のいずれか"))
            continue
        if st in NEEDS_SOURCE and not (under & NO_SOURCE_REQUIRED):
            if _empty(node.get("source")):
                issues.append(Issue(
                    "error", relfile, path or "status",
                    f"status: {st} なのに source が空",
                    "PMID:12345678 か DOI を入れる。出典が無いなら inferred に落とす"))
        if st == "figure_read":
            # どの図を読んだかが書かれていないと、後から検算できない。
            # これが無いと figure_read は「根拠を示さず verified を名乗る」抜け道になる。
            if _empty(node.get("figure")):
                issues.append(Issue(
                    "error", relfile, path or "status",
                    "status: figure_read なのに figure（どの図か）が空",
                    'その出典のどの図から読んだかを書く。例: figure: "Fig. 1"'))
            if not _empty(node.get("figure")) and _empty(node.get("read_note")):
                issues.append(Issue(
                    "warn", relfile, path or "status",
                    "図から何をどう読み取ったかの説明（read_note）が無い",
                    "目視解釈なので、読み取りの根拠を残しておくと検算できる"))
        if st == "inferred" and strict_inferred:
            issues.append(Issue("error", relfile, path or "status",
                                "strict_inferred が有効: inferred が残っている"))
        src = node.get("source")
        if isinstance(src, str) and src.strip():
            for one in [s.strip() for s in re.split(r"[;,]\s*", src) if s.strip()]:
                if not (PMID_RE.match(one) or DOI_RE.match(one) or one.startswith("http")):
                    issues.append(Issue(
                        "warn", relfile, f"{path}.source",
                        f"出典の書式が想定外: {one!r}",
                        "'PMID:26867490' / '10.1016/j.jhep.2015.10.002' / URL のいずれか。"
                        "タイトルだけの引用は不可（CLAUDE.md §0.1）"))


def collect_refs(obj, acc: set[str]):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == "source" and isinstance(v, str) and v.strip():
                for one in re.split(r"[;,]\s*", v):
                    if one.strip():
                        acc.add(one.strip())
            elif k == "refs" and isinstance(v, list):
                for r in v:
                    if isinstance(r, dict):
                        if r.get("pmid"):
                            acc.add(f"PMID:{r['pmid']}")
                        if r.get("doi"):
                            acc.add(str(r["doi"]))
            else:
                collect_refs(v, acc)
    elif isinstance(obj, list):
        for v in obj:
            collect_refs(v, acc)


# --------------------------------------------------------------------------
# 分子
# --------------------------------------------------------------------------
def check_molecule(mid: str, m: dict, relfile: str, ds: Dataset,
                   hgnc: dict, issues: list[Issue]):
    cells = set(ds.cells)
    fam_by_axis = {ax: {fid for fid, (f, _) in ds.families.items() if f.get("axis") == ax}
                   for ax in AXES}

    # --- gene / HGNC ---
    gene = m.get("gene")
    if _empty(gene):
        issues.append(Issue("error", relfile, "gene", "gene が無い（HGNC シンボルを入れる）"))
    elif hgnc:
        info = hgnc.get(str(gene))
        if info is None:
            issues.append(Issue("error", relfile, "gene",
                                f"'{gene}' は HGNC に無いシンボル"))
        elif info[0] != "approved":
            kind = "旧シンボル" if info[0] == "prev" else "別名(alias)"
            issues.append(Issue("error", relfile, "gene",
                                f"'{gene}' は{kind}であって承認シンボルではない",
                                f"承認シンボルは '{info[1]}'"))

    # --- 手書き禁止フィールド ---
    if "system_closure" in m:
        issues.append(Issue(
            "error", relfile, "system_closure",
            "system_closure は手で書かない（CLAUDE.md §5）",
            "build.py が producers/receivers から自動判定する。この行を削除する"))
    if "grp" in m or "family" in m:
        issues.append(Issue(
            "error", relfile, "family/grp",
            "単一の family / grp フィールドは使わない（v0 の設計ミス）",
            "families: {structure:, receptor:, naming:} の3軸で持つ"))

    # --- families 3軸 ---
    fams = m.get("families")
    if not isinstance(fams, dict):
        issues.append(Issue("error", relfile, "families",
                            "families が無い（structure / receptor / naming の3軸）"))
    else:
        for ax in AXES:
            if ax not in fams:
                issues.append(Issue("warn", relfile, f"families.{ax}",
                                    f"{ax} 軸が未設定"))
                continue
            fid = fams[ax]
            if _empty(fid):
                continue
            if fid not in ds.families:
                issues.append(Issue("error", relfile, f"families.{ax}",
                                    f"存在しない family_id '{fid}'"))
            elif fid not in fam_by_axis[ax]:
                actual = ds.families[fid][0].get("axis")
                issues.append(Issue("error", relfile, f"families.{ax}",
                                    f"'{fid}' は axis: {actual} のファミリー。{ax} 軸に置けない"))
            else:
                members = ds.families[fid][0].get("members") or []
                if gene and str(gene) not in members and mid not in members:
                    issues.append(Issue(
                        "warn", relfile, f"families.{ax}",
                        f"'{fid}' の members に {gene} が入っていない",
                        f"data/families/{ax}.yaml の members に追加すると両方向から辿れる"))

    # --- producers / receivers ---
    for side in ("producers", "receivers"):
        block = m.get(side)
        if block is None:
            issues.append(Issue("warn", relfile, side, f"{side} が無い"))
            continue
        if not isinstance(block, dict):
            issues.append(Issue("error", relfile, side,
                                f"{side} は {{cell_id: {{level, source, status}}}} 形式"))
            continue
        for cid, v in block.items():
            ppath = f"{side}.{cid}"
            if cid not in cells:
                issues.append(Issue("error", relfile, ppath,
                                    f"存在しない cell_id '{cid}'",
                                    "data/cells/ に登録するか、綴りを直す"))
            if not isinstance(v, dict):
                issues.append(Issue("error", relfile, ppath, "level/source/status のマッピングが必要"))
                continue
            lv = v.get("level")
            if not isinstance(lv, int) or isinstance(lv, bool) or not 0 <= lv <= 3:
                issues.append(Issue("error", relfile, f"{ppath}.level",
                                    f"level: {lv!r} は不正（0=なし 1=文脈依存 2=明確 3=主要）"))

    # --- effects ---
    eff = m.get("effects")
    if isinstance(eff, dict):
        for k in EFFECT_KEYS:
            if k not in eff:
                issues.append(Issue("warn", relfile, f"effects.{k}", f"effects.{k} が無い"))
                continue
            v = eff[k]
            if not isinstance(v, int) or isinstance(v, bool) or not -2 <= v <= 2:
                issues.append(Issue("error", relfile, f"effects.{k}",
                                    f"{v!r} は不正（-2〜+2 の整数）"))
        if isinstance(eff.get("note"), str) and "平均" in eff["note"]:
            pass
    elif eff is not None:
        issues.append(Issue("error", relfile, "effects", "effects はマッピング"))

    # --- form / range の語彙 ---
    form = m.get("form")
    if isinstance(form, dict) and not _empty(form.get("secretion")):
        if form["secretion"] not in SECRETION:
            issues.append(Issue("error", relfile, "form.secretion",
                                f"{form['secretion']!r} は不正",
                                " / ".join(sorted(SECRETION))))
    if isinstance(form, dict):
        if _empty(form.get("oligomer")):
            issues.append(Issue("warn", relfile, "form.oligomer",
                                "会合状態が未設定（図では単量体として描かれる）"))
        elif form["oligomer"] not in OLIGOMER:
            issues.append(Issue("error", relfile, "form.oligomer",
                                f"{form['oligomer']!r} は不正",
                                " / ".join(sorted(OLIGOMER))))
    rng = m.get("range")
    if isinstance(rng, dict) and not _empty(rng.get("mode")):
        if rng["mode"] not in RANGE_MODES:
            issues.append(Issue("error", relfile, "range.mode",
                                f"{rng['mode']!r} は不正",
                                " / ".join(sorted(RANGE_MODES))))

    # --- assays: プライマー ---
    assays = m.get("assays") or {}
    qpcr = assays.get("qpcr") or {}
    for i, pr in enumerate(qpcr.get("primers") or []):
        base = f"assays.qpcr.primers[{i}]"
        if not isinstance(pr, dict):
            issues.append(Issue("error", relfile, base, "プライマーはマッピング"))
            continue
        for key in ("fwd", "rev"):
            seq = pr.get(key)
            if _empty(seq):
                if pr.get("status") == "verified":
                    issues.append(Issue("error", relfile, f"{base}.{key}",
                                        f"verified なのに {key} が空"))
                continue
            s = str(seq).strip().upper()
            shown = {" ": "空白", "\t": "タブ", "\n": "改行", "-": "ハイフン"}
            bad = sorted(set(re.sub(r"[ATGC]", "", s)))
            if bad:
                issues.append(Issue("error", relfile, f"{base}.{key}",
                                    "配列に ATGC 以外の文字: "
                                    + " ".join(shown.get(ch, repr(ch)) for ch in bad),
                                    "縮重塩基やスペースも不可。一次情報から貼り直す"))
            elif not 15 <= len(s) <= 35:
                issues.append(Issue("warn", relfile, f"{base}.{key}",
                                    f"長さ {len(s)} nt は qPCR プライマーとして異例",
                                    "貼り間違いでないか確認する"))
        if pr.get("status") == "verified" and _empty(pr.get("source")):
            issues.append(Issue("error", relfile, f"{base}.source",
                                "verified なのに出典が無い（§0.1 プライマー配列は推測禁止）"))

    for i, ab in enumerate(assays.get("antibody") or []):
        if isinstance(ab, dict) and ab.get("status") == "verified":
            if _empty(ab.get("clone")) and _empty(ab.get("cat")):
                issues.append(Issue("error", relfile, f"assays.antibody[{i}]",
                                    "verified なのにクローン番号もカタログ番号も無い"))

    # --- own_data ---
    for i, od in enumerate(m.get("own_data") or []):
        base = f"own_data[{i}]"
        if not isinstance(od, dict):
            issues.append(Issue("error", relfile, base, "own_data の要素はマッピング"))
            continue
        if _empty(od.get("date")):
            issues.append(Issue("error", relfile, f"{base}.date",
                                "own_data には日付が要る（いつの実測か追えなくなる）"))
        elif not re.match(r"^\d{4}-\d{2}-\d{2}$", str(od["date"])):
            issues.append(Issue("warn", relfile, f"{base}.date",
                                f"{od['date']!r} は YYYY-MM-DD 形式でない"))
        if od.get("value") is not None and _empty(od.get("unit")):
            issues.append(Issue("error", relfile, f"{base}.unit",
                                "value があるのに unit が無い"))

    # --- 分子間辺 ---
    for i, rel in enumerate(m.get("relations") or []):
        base = f"relations[{i}]"
        if not isinstance(rel, dict):
            issues.append(Issue("error", relfile, base, "relations の要素はマッピング"))
            continue
        rt = rel.get("type")
        if rt not in RELATION_TYPES:
            issues.append(Issue("error", relfile, f"{base}.type",
                                f"{rt!r} は不正",
                                " / ".join(sorted(RELATION_TYPES))))
        if rt == "shares_receptor_with":
            issues.append(Issue(
                "error", relfile, f"{base}.type",
                "受容体の共有は手で書かない",
                "受容体軸のファミリー所属から build.py が導出する"))
        tgt = rel.get("target")
        if _empty(tgt):
            issues.append(Issue("error", relfile, f"{base}.target", "target が無い"))
        elif tgt not in ds.molecules:
            issues.append(Issue("error", relfile, f"{base}.target",
                                f"存在しない分子 '{tgt}'"))
        elif tgt == mid:
            issues.append(Issue("error", relfile, f"{base}.target",
                                "自分自身に関係を張っている"))

    # --- pathway / theme 参照 ---
    if "pathways" in m:
        issues.append(Issue(
            "error", relfile, "pathways",
            "分子側に pathways を手で書かない",
            "data/pathways/ の steps[].molecules から build.py が導出する"))
    for i, ref in enumerate(m.get("themes") or []):
        rid = ref.get("id") if isinstance(ref, dict) else ref
        if rid and rid not in ds.themes:
            issues.append(Issue("error", relfile, f"themes[{i}]",
                                f"存在しない theme_id '{rid}'"))


# --------------------------------------------------------------------------
# 細胞 / ファミリー
# --------------------------------------------------------------------------
def check_cells(ds: Dataset, issues: list[Issue]):
    for cid, (c, relfile) in ds.cells.items():
        if _empty(c.get("label_ja")):
            issues.append(Issue("warn", relfile, f"{cid}.label_ja", "label_ja が無い"))
        if "kind" in c:
            issues.append(Issue(
                "error", relfile, f"{cid}.kind",
                "kind は使わない（分化段階と由来という無関係な2軸を1フィールドに潰していた）",
                "stage: pluripotent|progenitor|mature と origin: primary|line|ipsc_derived に分ける"))
        stage = c.get("stage")
        if stage is None:
            issues.append(Issue("error", relfile, f"{cid}.stage", "stage が無い",
                                " / ".join(sorted(CELL_STAGES))))
        elif stage not in CELL_STAGES:
            issues.append(Issue("error", relfile, f"{cid}.stage",
                                f"{stage!r} は不正", " / ".join(sorted(CELL_STAGES))))
        if "origin" not in c:
            issues.append(Issue(
                "error", relfile, f"{cid}.origin", "origin が無い",
                "primary / line / ipsc_derived。培養していない参照ノードなら null を明示する"))
        elif c["origin"] is not None and c["origin"] not in CELL_ORIGINS:
            issues.append(Issue("error", relfile, f"{cid}.origin",
                                f"{c['origin']!r} は不正",
                                " / ".join(sorted(CELL_ORIGINS)) + " / null"))
        # 系内の細胞は実際に培養しているはずなので origin が要る
        if c.get("in_system", True) and c.get("origin") is None and "origin" in c:
            issues.append(Issue("error", relfile, f"{cid}.origin",
                                "系内の細胞なのに origin が null",
                                "系内＝実際に培養している細胞。由来を書く"))
        if not _empty(c.get("label_generic")) and _empty(c.get("label_short")):
            issues.append(Issue("warn", relfile, f"{cid}.label_short",
                                "総称はあるが実体の呼称が無い"))
        color = c.get("color")
        if _empty(color):
            issues.append(Issue("warn", relfile, f"{cid}.color",
                                "color が無い（図の色が一貫しなくなる）"))
        elif not re.match(r"^#[0-9A-Fa-f]{6}$", str(color)):
            issues.append(Issue("error", relfile, f"{cid}.color",
                                f"{color!r} は #RRGGBB 形式でない"))
        for key in ("source_of",):
            v = c.get(key)
            if v and v not in ds.cells:
                issues.append(Issue("error", relfile, f"{cid}.{key}",
                                    f"存在しない cell_id '{v}'"))
        diff = c.get("differentiation")
        if isinstance(diff, dict):
            frm = diff.get("from")
            if frm and frm not in ds.cells:
                issues.append(Issue("error", relfile, f"{cid}.differentiation.from",
                                    f"存在しない cell_id '{frm}'"))
            if frm == cid:
                issues.append(Issue("error", relfile, f"{cid}.differentiation.from",
                                    "自分自身から分化することになっている"))
            for mk in diff.get("markers_gained") or []:
                pass  # マーカーは遺伝子とは限らない（CD 番号等）ので HGNC 照合しない

    # 分化グラフの循環
    parent = {cid: (c.get("differentiation") or {}).get("from")
              for cid, (c, _) in ds.cells.items()}
    for cid in ds.cells:
        seen, cur = [], cid
        while cur:
            if cur in seen:
                relfile = ds.cells[cid][1]
                issues.append(Issue("error", relfile, f"{cid}.differentiation",
                                    f"分化グラフが循環している: {' → '.join(seen + [cur])}"))
                break
            seen.append(cur)
            cur = parent.get(cur)


def check_families(ds: Dataset, issues: list[Issue]):
    genes = {str(m.get("gene")) for m, _ in ds.molecules.values() if m.get("gene")}
    ids = set(ds.molecules)
    for fid, (f, relfile) in ds.families.items():
        ax = f.get("axis")
        if ax not in AXES:
            issues.append(Issue("error", relfile, f"{fid}.axis",
                                f"{ax!r} は不正", " / ".join(sorted(AXES))))
        par = f.get("parent")
        if par:
            if par not in ds.families:
                issues.append(Issue("error", relfile, f"{fid}.parent",
                                    f"存在しない family_id '{par}'"))
            elif ds.families[par][0].get("axis") != ax:
                issues.append(Issue("error", relfile, f"{fid}.parent",
                                    f"親 '{par}' は別の軸（{ds.families[par][0].get('axis')}）"))
        for i, ch in enumerate(f.get("children") or []):
            if ch not in ds.families:
                issues.append(Issue("error", relfile, f"{fid}.children[{i}]",
                                    f"存在しない family_id '{ch}'"))
            elif ds.families[ch][0].get("parent") != fid:
                issues.append(Issue("warn", relfile, f"{fid}.children[{i}]",
                                    f"'{ch}' の parent が '{fid}' を指していない"))
        for i, mem in enumerate(f.get("mismatch") or []):
            if isinstance(mem, dict) and mem.get("member"):
                if mem["member"] not in genes and mem["member"] not in ids:
                    issues.append(Issue(
                        "warn", relfile, f"{fid}.mismatch[{i}].member",
                        f"'{mem['member']}' の分子ファイルがまだ無い",
                        "data/molecules/ に追加すると相互リンクが張れる"))
        # members に未登録の分子が並ぶのは正常（これから書く分）。WARN も出さない。

    # 親子の循環
    for fid in ds.families:
        seen, cur = [], fid
        while cur:
            if cur in seen:
                issues.append(Issue("error", ds.families[fid][1], f"{fid}.parent",
                                    f"ファミリーの親子が循環: {' → '.join(seen + [cur])}"))
                break
            seen.append(cur)
            cur = ds.families[cur][0].get("parent") if cur in ds.families else None


def check_architectures(ds: Dataset, issues: list[Issue]):
    for aid, (a, relfile) in ds.architectures.items():
        if a.get("tm") not in TM_COUNTS:
            issues.append(Issue("error", relfile, f"{aid}.tm",
                                f"{a.get('tm')!r} は不正（膜貫通回数）",
                                " / ".join(str(x) for x in sorted(TM_COUNTS, key=str))))
        ch = a.get("chains")
        if not isinstance(ch, int) or isinstance(ch, bool) or ch < 1:
            issues.append(Issue("error", relfile, f"{aid}.chains",
                                f"{ch!r} は不正（1以上の整数）"))
        if _empty(a.get("crossing")):
            issues.append(Issue("error", relfile, f"{aid}.crossing",
                                "膜を越えてシグナルを伝える手段が書かれていない",
                                "形だけ書いても『なぜその形なのか』が伝わらない"))
        if a.get("tm") == 1 and a.get("chains") == 1:
            issues.append(Issue(
                "warn", relfile, f"{aid}",
                "1回膜貫通で1本鎖だと、単独では細胞内へ構造変化を伝えられない",
                "二量体化するなら chains を2以上にする。例外なら crossing にその旨を書く"))


def check_production(ds: Dataset, issues: list[Issue]):
    """産生経路。受容側とは別の語彙で検証する。"""
    seen_mol: dict[str, str] = {}
    for pid, (pr, relfile) in ds.production.items():
        mol = pr.get("molecule")
        if _empty(mol):
            issues.append(Issue("error", relfile, f"{pid}.molecule",
                                "どの分子の産生経路かが無い"))
        elif mol not in ds.molecules:
            issues.append(Issue("error", relfile, f"{pid}.molecule",
                                f"存在しない分子 '{mol}'"))
        elif mol in seen_mol:
            issues.append(Issue("warn", relfile, f"{pid}.molecule",
                                f"'{mol}' の産生経路が複数ある（既出: {seen_mol[mol]}）"))
        else:
            seen_mol[str(mol)] = pid
        if _empty(pr.get("outcome")):
            issues.append(Issue("warn", relfile, f"{pid}.outcome",
                                "この経路の要点（何が律速か）が書かれていない"))

        steps = pr.get("steps")
        if not isinstance(steps, list) or not steps:
            issues.append(Issue("error", relfile, f"{pid}.steps", "steps が無い"))
            continue
        kinds, signals, ids = [], [], set()
        for i, st in enumerate(steps):
            base = f"{pid}.steps[{i}]"
            if not isinstance(st, dict):
                issues.append(Issue("error", relfile, base, "ステップはマッピング"))
                continue
            sid = st.get("id")
            if _empty(sid):
                issues.append(Issue("error", relfile, f"{base}.id", "id が無い"))
            elif sid in ids:
                issues.append(Issue("error", relfile, f"{base}.id",
                                    f"ステップ id '{sid}' が重複している"))
            else:
                ids.add(sid)
            k = st.get("kind")
            kinds.append(k)
            if k not in PROD_STEP_KINDS:
                issues.append(Issue("error", relfile, f"{base}.kind",
                                    f"{k!r} は不正",
                                    " / ".join(sorted(PROD_STEP_KINDS))))
            if _empty(st.get("label")):
                issues.append(Issue("error", relfile, f"{base}.label", "label が無い"))
            if k == "trigger":
                sg = st.get("signal")
                if not isinstance(sg, int) or isinstance(sg, bool) or sg < 1:
                    issues.append(Issue(
                        "error", relfile, f"{base}.signal",
                        "trigger には何番目のシグナルかを書く（1以上の整数）",
                        "2シグナル制御かどうかはこの番号から導出する"))
                else:
                    signals.append(sg)
            elif st.get("signal") is not None:
                issues.append(Issue("error", relfile, f"{base}.signal",
                                    "signal は trigger のステップにだけ付ける"))
            rt = st.get("route")
            if rt is not None and rt not in RELEASE_ROUTES:
                issues.append(Issue("error", relfile, f"{base}.route",
                                    f"{rt!r} は不正", " / ".join(sorted(RELEASE_ROUTES))))
            if rt and k != "release":
                issues.append(Issue("error", relfile, f"{base}.route",
                                    "route は release のステップにだけ付ける"))
        if "release" not in kinds:
            issues.append(Issue("error", relfile, f"{pid}.steps",
                                "release のステップが無い（どうやって外に出るかが不明）"))
        if "transcription" not in kinds:
            issues.append(Issue("warn", relfile, f"{pid}.steps",
                                "transcription のステップが無い"))


def check_pathways(ds: Dataset, issues: list[Issue]):
    fam_receptor = {fid for fid, (f, _) in ds.families.items()
                    if f.get("axis") == "receptor"}
    for pid, (pw, relfile) in ds.pathways.items():
        fam = pw.get("family")
        if _empty(fam):
            issues.append(Issue("warn", relfile, f"{pid}.family",
                                "対応する受容体ファミリーが未設定"))
        elif fam not in ds.families:
            issues.append(Issue("error", relfile, f"{pid}.family",
                                f"存在しない family_id '{fam}'"))
        elif fam not in fam_receptor:
            issues.append(Issue("error", relfile, f"{pid}.family",
                                f"'{fam}' は受容体軸のファミリーではない",
                                "経路は受容体軸のファミリーに対応させる"))
        if _empty(pw.get("outcome")):
            issues.append(Issue("warn", relfile, f"{pid}.outcome",
                                "この経路が何をするかが書かれていない"))

        steps = pw.get("steps")
        if not isinstance(steps, list) or not steps:
            issues.append(Issue("error", relfile, f"{pid}.steps", "steps が無い"))
            continue
        seen_step, seen_mol = set(), {}
        for i, st in enumerate(steps):
            base = f"{pid}.steps[{i}]"
            if not isinstance(st, dict):
                issues.append(Issue("error", relfile, base, "ステップはマッピング"))
                continue
            sid = st.get("id")
            if _empty(sid):
                issues.append(Issue("error", relfile, f"{base}.id", "ステップに id が無い"))
            elif sid in seen_step:
                issues.append(Issue("error", relfile, f"{base}.id",
                                    f"ステップ id '{sid}' が重複している"))
            else:
                seen_step.add(sid)
            arch = st.get("architecture")
            if arch is not None and arch not in ds.architectures:
                issues.append(Issue("error", relfile, f"{base}.architecture",
                                    f"存在しない architecture '{arch}'",
                                    "data/architectures/registry.yaml に登録する"))
            if arch and st.get("kind") != "receptor":
                issues.append(Issue("error", relfile, f"{base}.architecture",
                                    "architecture は receptor のステップにだけ付ける"))
            if st.get("kind") == "receptor" and not arch:
                issues.append(Issue("warn", relfile, f"{base}.architecture",
                                    "膜での形が未設定（図では一般的な2本鎖として描かれる）"))
            if st.get("kind") not in STEP_KINDS:
                issues.append(Issue("error", relfile, f"{base}.kind",
                                    f"{st.get('kind')!r} は不正",
                                    " / ".join(sorted(STEP_KINDS))))
            if _empty(st.get("label")):
                issues.append(Issue("error", relfile, f"{base}.label", "label が無い"))
            for field in ("molecules", "induces"):
                for j, sym in enumerate(st.get(field) or []):
                    if sym not in ds.molecules:
                        issues.append(Issue("error", relfile, f"{base}.{field}[{j}]",
                                            f"存在しない分子 '{sym}'"))
                    elif field == "molecules":
                        if sym in seen_mol:
                            issues.append(Issue(
                                "error", relfile, f"{base}.molecules[{j}]",
                                f"'{sym}' が同じ経路の複数ステップに現れる"
                                f"（既出: {seen_mol[sym]}）"))
                        else:
                            seen_mol[sym] = sid
        kinds = [st.get("kind") for st in steps if isinstance(st, dict)]
        if "ligand" not in kinds:
            issues.append(Issue("warn", relfile, f"{pid}.steps",
                                "ligand のステップが無い（入口が不明）"))
        if "output" not in kinds:
            issues.append(Issue("warn", relfile, f"{pid}.steps",
                                "output のステップが無い（この経路が何を出すか不明）"))


def check_config(ds: Dataset, issues: list[Issue]):
    sysc = (ds.config.get("system") or {}).get("in_system") or []
    for cid in sysc:
        if cid not in ds.cells:
            issues.append(Issue("error", "config.yaml", "system.in_system",
                                f"存在しない cell_id '{cid}'"))
    mat = (ds.config.get("matrix") or {}).get("cells") or []
    for cid in mat:
        if cid not in ds.cells:
            issues.append(Issue("error", "config.yaml", "matrix.cells",
                                f"存在しない cell_id '{cid}'"))
    if mat and len(mat) != len(set(mat)):
        issues.append(Issue("error", "config.yaml", "matrix.cells", "細胞が重複している"))
    if "include_own_data" not in ds.config:
        issues.append(Issue("warn", "config.yaml", "include_own_data",
                            "include_own_data が無い（安全側の false として扱う）"))
    # cells 側の in_system と config の一覧が食い違っていないか
    for cid, (c, relfile) in ds.cells.items():
        if "in_system" in c and sysc:
            declared = bool(c["in_system"])
            listed = cid in sysc
            if declared != listed:
                issues.append(Issue(
                    "error", relfile, f"{cid}.in_system",
                    f"cells は in_system: {declared} だが config.yaml の一覧では {listed}",
                    "どちらかに揃える。食い違うと system_closure の判定が壊れる"))


# --------------------------------------------------------------------------
# 出典の実在確認（--check-refs）
# --------------------------------------------------------------------------
def check_refs_online(refs: set[str], issues: list[Issue], where: dict[str, list[str]]):
    import urllib.parse
    import urllib.request

    pmids = sorted({m.group(1) for r in refs if (m := PMID_RE.match(r))})
    if pmids:
        found: dict[str, str] = {}
        for i in range(0, len(pmids), 100):
            chunk = pmids[i:i + 100]
            url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                   f"?db=pubmed&retmode=json&id={','.join(chunk)}")
            try:
                with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310
                    res = json.loads(r.read().decode())["result"]
            except Exception as e:  # pragma: no cover
                issues.append(Issue("warn", "-", "--check-refs",
                                    f"PubMed に問い合わせできなかった: {e}"))
                return
            for uid in res.get("uids", []):
                rec = res[uid]
                # 存在しない uid でも esummary は殻のレコードを返す。error 付きは「無い」扱い。
                if rec.get("error") or not rec.get("title"):
                    continue
                found[uid] = (f"{rec.get('sortfirstauthor')} / {rec.get('source')} "
                              f"{rec.get('pubdate')} / {rec.get('title', '')[:70]}")
        for pmid in pmids:
            key = f"PMID:{pmid}"
            locs = ", ".join(sorted(set(where.get(key, [])))[:3])
            if pmid not in found:
                issues.append(Issue("error", locs or "-", key,
                                    "この PMID は PubMed に存在しない",
                                    "番号の写し間違いか、実在しない文献（CLAUDE.md §0.1）"))
            else:
                issues.append(Issue("warn", locs or "-", key,
                                    f"実在を確認: {found[pmid]}",
                                    "内容が主張と一致しているか目視で確かめる"))

    dois = sorted({r for r in refs if DOI_RE.match(r)})
    for doi in dois:
        clean = re.sub(r"^DOI:\s*", "", doi, flags=re.I)
        url = "https://api.crossref.org/works/" + urllib.parse.quote(clean, safe="")
        locs = ", ".join(sorted(set(where.get(doi, [])))[:3])
        try:
            with urllib.request.urlopen(url, timeout=30) as r:  # noqa: S310
                msg = json.loads(r.read().decode())["message"]
            title = (msg.get("title") or [""])[0][:70]
            issues.append(Issue("warn", locs or "-", doi, f"実在を確認: {title}"))
        except Exception:
            issues.append(Issue("error", locs or "-", doi, "この DOI は Crossref で解決できない"))


# --------------------------------------------------------------------------
def run(strict: bool, check_refs: bool) -> tuple[list[Issue], Dataset]:
    issues: list[Issue] = []
    ds = load_dataset(issues)
    hgnc = load_hgnc(ds, issues)
    strict_inferred = bool((ds.config.get("validate") or {}).get("strict_inferred"))

    check_config(ds, issues)
    check_cells(ds, issues)
    check_families(ds, issues)
    check_architectures(ds, issues)
    check_pathways(ds, issues)
    check_production(ds, issues)

    all_refs: set[str] = set()
    ref_where: dict[str, list[str]] = {}
    for store in (ds.molecules, ds.cells, ds.families, ds.pathways, ds.themes,
                  ds.architectures, ds.production):
        for oid, (obj, relfile) in store.items():
            check_status_and_sources(obj, relfile, issues, strict_inferred)
            acc: set[str] = set()
            collect_refs(obj, acc)
            for r in acc:
                ref_where.setdefault(r, []).append(relfile)
            all_refs |= acc

    for mid, (m, relfile) in ds.molecules.items():
        check_molecule(mid, m, relfile, ds, hgnc, issues)

    if check_refs:
        check_refs_online(all_refs, issues, ref_where)

    if strict:
        for i in issues:
            if i.level == "warn":
                i.level = "error"
    return issues, ds


def main() -> int:
    ap = argparse.ArgumentParser(description="data/ の整合性チェック")
    ap.add_argument("--strict", action="store_true", help="WARN もエラー扱いにする")
    ap.add_argument("--check-refs", action="store_true",
                    help="PMID/DOI が実在するか問い合わせる（要ネットワーク）")
    ap.add_argument("--json", action="store_true", help="JSON で出力")
    args = ap.parse_args()

    issues, ds = run(args.strict, args.check_refs)
    errors = [i for i in issues if i.level == "error"]
    warns = [i for i in issues if i.level == "warn"]

    if args.json:
        print(json.dumps({
            "ok": not errors,
            "counts": {"error": len(errors), "warn": len(warns)},
            "issues": [i.__dict__ for i in issues],
        }, ensure_ascii=False, indent=2))
        return 1 if errors else 0

    print(f"分子 {len(ds.molecules)} / 細胞 {len(ds.cells)} / ファミリー {len(ds.families)} "
          f"/ 受容経路 {len(ds.pathways)} / 産生経路 {len(ds.production)} "
          f"/ 受容体の形 {len(ds.architectures)} "
          f"/ テーマ {len(ds.themes)}")
    if errors:
        print(f"\n■ ERROR {len(errors)} 件（ビルドを止める）")
        for i in errors:
            print(i.fmt())
    if warns:
        print(f"\n■ WARN {len(warns)} 件")
        for i in warns:
            print(i.fmt())
    if not issues:
        print("\n問題なし。")
    print()
    if errors:
        print(f"検証失敗: ERROR {len(errors)} 件。直すまでビルドしない。")
        return 1
    print(f"検証通過（WARN {len(warns)} 件）。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
