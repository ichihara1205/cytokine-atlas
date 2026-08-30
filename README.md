# Cytokine Atlas

肝共培養系（および iPS 由来細胞の分化・成熟）を対象にしたサイトカイン参照サイト。
`data/` の YAML からサイトを生成する。**資産はデータであってページではない。**

## 使い方

```bash
python3 -m pip install pyyaml
python3 scripts/validate.py          # data/ の整合性チェック（ERROR があれば exit 1）
python3 scripts/build.py --serve     # site/ を生成してローカルサーバを立てる
```

補助:

```bash
python3 scripts/validate.py --check-refs   # PMID/DOI が実在するか PubMed/Crossref に問い合わせ
python3 scripts/validate.py --strict       # WARN もエラー扱い（公開前の総点検）
python3 scripts/fetch_hgnc.py              # HGNC の遺伝子シンボル照合表を更新（半年に一度程度）
```

## 記述の信頼性

**このサイトの記述の大半は未検証。**

| status | 意味 |
|---|---|
| `verified` | 一次文献または本人の実測で確認済み。出典あり |
| `inferred` | 一般知識から書いたもの。**未検証** |
| `todo` | 未着手 |

`validate.py` が以下を機械的に保証する（違反はビルド失敗）:

- `verified` なのに `source` が空
- `gene` が HGNC 承認シンボルでない（公式データ 45,031 件と照合）
- プライマー配列に ATGC 以外の文字
- 存在しない `cell_id` / `family_id` / 分子 への参照
- 導出すべきフィールドの手書き（`system_closure` / 分子側の `pathways` / 受容体の共有）

残作業は生成サイトの `/todo` に全件出る。

## 構成

| ディレクトリ | 中身 |
|---|---|
| `data/molecules/` | 分子（1分子1ファイル） |
| `data/cells/` | 細胞レジストリ。`stage`（分化段階）と `origin`（由来）は別軸 |
| `data/families/` | 3軸のファミリー（構造相同性 / 受容体共有 / 慣用名） |
| `data/pathways/` | 受容側の経路（外 → 中） |
| `data/production/` | 産生側の経路（中 → 外） |
| `data/architectures/` | 受容体の膜での形と、その形が要求する伝達手段 |
| `scripts/` | `validate.py` / `build.py` / `fetch_hgnc.py` |
| `site/` | 生成物。**手で編集しない**（gitignore 済み） |

図はすべてデータから生成した inline SVG。外部の図は一切使っていない。

詳しい規約は `CLAUDE.md`、変更履歴は `CHANGELOG.md`。

## 公開について

GitHub Actions が `main` への push でビルドして GitHub Pages に配信する。
**リポジトリが private でも、公開されたサイトは URL を知っていれば誰でも見られる。**
`own_data`（未発表の実測）は `config.yaml` の `include_own_data: false` により
HTML にも JSON にも出力されない。
