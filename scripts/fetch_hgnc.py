#!/usr/bin/env python3
"""HGNC の公式全遺伝子セットを取得し、validate.py が使う軽量な照合表を作る。

出力: data/hgnc/symbols.tsv   （symbol \t kind \t current_symbol）
  kind = approved | prev | alias

approved 以外を持たせているのは、validate.py が
「IL8 は旧シンボル。現行は CXCL8」と直せる形でエラーを出すため。

使い方:
    python3 scripts/fetch_hgnc.py            # ダウンロードして再生成
    python3 scripts/fetch_hgnc.py --from FILE # 手元の hgnc_complete_set.txt から生成

HGNC のシンボルは更新される。半年に一度くらい流し直すこと。
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "hgnc" / "symbols.tsv"

URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


def load(src: str | None) -> str:
    if src:
        return Path(src).read_text(encoding="utf-8")
    sys.stderr.write(f"downloading {URL}\n")
    with urllib.request.urlopen(URL, timeout=300) as r:  # noqa: S310
        return r.read().decode("utf-8")


def build(text: str) -> list[tuple[str, str, str]]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    rows: list[tuple[str, str, str]] = []
    approved: set[str] = set()
    secondary: dict[str, tuple[str, str]] = {}

    for rec in reader:
        if rec.get("status") != "Approved":
            continue
        sym = (rec.get("symbol") or "").strip()
        if not sym:
            continue
        approved.add(sym)
        for field, kind in (("prev_symbol", "prev"), ("alias_symbol", "alias")):
            raw = (rec.get(field) or "").strip().strip('"')
            for other in filter(None, (s.strip() for s in raw.split("|"))):
                # prev を alias より優先。先勝ちにはしない
                if other not in secondary or kind == "prev":
                    secondary[other] = (kind, sym)

    for sym in sorted(approved):
        rows.append((sym, "approved", sym))
    for sym, (kind, current) in sorted(secondary.items()):
        if sym in approved:
            continue  # 承認済みシンボルが別遺伝子の alias でもある場合は承認を優先
        rows.append((sym, kind, current))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", help="ローカルの hgnc_complete_set.txt")
    args = ap.parse_args()

    rows = build(load(args.src))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["symbol", "kind", "current_symbol"])
        w.writerows(rows)

    n_app = sum(1 for r in rows if r[1] == "approved")
    print(f"{OUT.relative_to(ROOT)}: {len(rows)} 行 (approved {n_app})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
