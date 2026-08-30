# Cytokine Atlas — プロジェクト仕様

肝共培養系（および iPS 由来細胞の分化・成熟）を対象にしたサイトカイン参照サイト。
データファイル群から静的サイトを生成する。**資産はデータであってページではない。**

このファイルは Claude Code が毎セッション最初に読む。以下のルールは例外なく適用する。

---

## 0. Claude への厳格ルール（最重要）

### 0.1 捏造の禁止
- **プライマー配列、抗体クローン番号、カタログ番号、ロット番号を推測で書いてはならない。**
  出典が特定できない場合は該当フィールドを空にし、`status: todo` を立てる。
  「たぶんこれ」で埋めることは、このプロジェクトにおける最大の失敗。
- 数値（濃度、Cq、pg/mL など）も同様。ユーザーの実測か、原著の記載のみ。
- 文献は DOI か PMID を必ず添える。タイトルだけの引用は不可。存在を確認できない文献は書かない。

### 0.2 status を必ず付ける
各ブロックに以下のいずれかを付ける。省略不可。

| status | 意味 | サイト上の表示 |
|---|---|---|
| `verified` | 一次文献または本人の実測で確認済み。出典あり | 通常表示 |
| `inferred` | Claude が一般知識から書いた。未検証 | グレー＋点線枠＋「未検証」バッジ |
| `todo` | 未着手・要調査 | 空欄＋TODOリストに集計 |

- 新規に書いたものは原則 `inferred` から始める。`verified` に上げるのはユーザーの指示があったときだけ。
- `assays.qpcr.primers` と `assays.antibody` は **`verified` 以外はサイトに出力しない**（ビルドで落とす）。

### 0.3 ビルド前チェック（scripts/validate.py）
以下に該当したらビルドを失敗させる。
- `status: verified` なのに `source` が空
- 存在しない `cell_id` / `pathway_id` / `theme_id` を参照している
- `gene` が HGNC シンボルと一致しない
- プライマー配列に ATGC 以外の文字が入っている

### 0.4 変更の記録
データを変更したら `CHANGELOG.md` に1行足す（日付・分子・何を変えたか・出典）。
「いつ誰が何を根拠に書いたか」が追えなくなった時点でこのサイトは信用を失う。

---

## 1. フェーズ（この順で作る）

| # | 内容 | 完了条件 |
|---|---|---|
| **1** | **産生↔受容マップ ＋ ファミリー体系** | 分子50、細胞レジストリ確立、3軸ファミリー、5×5行列と絞り込みが動く |
| 2 | 経路図（ファミリー単位） | 主要8ファミリーの SVG 図。分子ページから相互リンク |
| 3 | テーマ別の全体像 | 肝の炎症 / 脂質代謝 / iKC成熟 / iPS→肝分化 の4枚 |
| 4 | 検証法DB | qPCR プライマー・抗体・RNA-seq。**verified のみ** |

フェーズ1が終わるまで、2以降のファイルは作らない。中途半端に4つ並走させると全部が未完成のまま重くなる。

---

## 2. 公開判断

当面は**自分専用**（ローカル or private repo）。公開はフェーズ1完了時点で判断する。

公開に切り替えるとき必要になるのは以下だけ。この3つを常に満たしておけば、判断は後回しにできる。
1. すべての `verified` に出典がある
2. 他社の図（R&D Systems など）を一切コピーしていない。**図はすべて自作 SVG**
3. `own_data`（本人の未発表実測）を出力から除外するフラグが効く → `config.yaml` の `include_own_data: false`

---

## 3. ディレクトリ構成

```
cytokine-atlas/
├── CLAUDE.md              # このファイル
├── CHANGELOG.md
├── config.yaml            # サイト設定・公開フラグ
├── data/
│   ├── molecules/         # IL6.yaml, TGFB1.yaml, ...  （1分子1ファイル）
│   ├── cells/             # 細胞レジストリ（後述）
│   ├── pathways/          # フェーズ2
│   └── themes/            # フェーズ3
├── scripts/
│   ├── validate.py        # 0.3 のチェック
│   └── build.py           # data/ → site/ を生成
├── site/                  # 生成物。手で編集しない
└── sources/               # PDF・スクショ・出典メモのキャッシュ
```

**site/ を直接編集してはならない。** 直したいものがあれば data/ か scripts/ を直す。

---

## 4. 細胞レジストリ（data/cells/）

細胞は5種固定ではない。**分化段階もノードとして持つ。**

```yaml
id: ikc
label_ja: iPS由来クッパー細胞
label_short: iKC
kind: mature            # mature | progenitor | line
lineage: myeloid
color: "#4FE08C"        # 図で一貫して使う色
source_of: ipsc         # 分化元（分化グラフの辺）
differentiation:
  from: premac
  protocol: "Tasnim et al. 2019 stage 3"
  markers_gained: [CD163, VSIG4, ID3]
  source: "PMID:30690240"
  status: verified
notes: "..."
```

### 辺の種類を混ぜないこと
- **通信辺**（molecule の producers/receivers から生成）＝ 成熟細胞どうしの空間的な会話
- **分化辺**（cells の `differentiation.from`）＝ 同じ細胞の時間軸

この2つは別のグラフとして描画する。同じ図に混ぜると意味が壊れる。

初期登録：`hepatocyte_pxb`, `hsc_lx2`, `hsc_qhsc`, `lsec_tmnk1`, `lsec_ilsec`, `kc_ikc`, `kc_pkc`, `mac_thp1`, `premac`, `ipsc`, `hepatoblast`, `neutrophil`(系外), `t_nk`(系外), `adipocyte`(系外)

「系外」の細胞も登録する。**受け手が系内にいないことを明示するのがこのサイトの中心的な価値**なので、系外ノードがないと判定できない。

---

## 4b. ファミリー体系（data/families/）

サイトカインの命名には**3つの互いに無関係な論理**が混在している。1本の木では表せない。
**必ず3軸を別々に持ち、UI で切り替える。単一の `family` フィールドに潰してはならない。**

| 軸 | id | 何を表すか | 代表例 |
|---|---|---|---|
| 構造相同性 | `structure` | 配列・フォールドの近さ | TGF-βスーパーファミリー（TGF-β/BMP/GDF/アクチビン/AMH/nodal/ミオスタチン を包含） |
| 受容体共有 | `receptor` | 受容体鎖・下流の共有 | gp130共有群（IL-6/IL-11/OSM/LIF/CNTF/CT-1）、γc群、IL10RB共有群 |
| 慣用名・通し番号 | `naming` | 人が検索で最初に叩く軸 | IL番号、CCL/CXCL番号、FGF番号 |

### 必ず表示すべき「名前と実体のズレ」
これを可視化することがこのレイヤの存在意義。`mismatch` フィールドに明示する。

- **IL-28A / IL-28B / IL-29 = IFN-λ2 / λ3 / λ1**。名前は IL、機能は III型 IFN、受容体は IFNLR1 + **IL10RB**（IL-10ファミリーと鎖を共有）
- **BMP・GDF・アクチビン・ミオスタチン(GDF8)・AMH・nodal は全部 TGF-βスーパーファミリー**。SMAD1/5/8 枝と SMAD2/3 枝に分岐するだけ
- **GDF15 は構造上 TGF-β系だが受容体は GFRAL+RET**。SMAD を使わない完全な例外
- **IL の番号は発見順**であり系統を反映しない。IL-1/18/33/36/37/38 が一族、IL-2/4/7/9/15/21 が γc 一族
- **OSM/LIF/CNTF/CT-1/IL-11 は IL-6 と同族**（gp130 共有）だが名前が揃っていない
- **FGF19/21/23 だけが内分泌型 FGF**。ヘパリン結合能を失い α/β-Klotho 依存になる

### スキーマ
```yaml
id: tgf-beta-superfamily
axis: structure            # structure | receptor | naming
label_ja: TGF-βスーパーファミリー
parent: null
children: [tgfb-smad23-branch, bmp-smad158-branch]
members: [TGFB1, TGFB2, TGFB3, BMP2, BMP4, GDF2, GDF8, GDF15, INHBA, AMH, NODAL]
defining_feature: "システインノット構造。二量体で I型/II型セリンスレオニンキナーゼ受容体に結合"
mismatch:
  - member: GDF15
    note: "構造は本ファミリーだが受容体は GFRAL+RET。SMAD 非依存"
    source: "PMID:28965740"
    status: verified
source: "PMID:27141051"
status: verified
```

分子側は `families: {structure: ..., receptor: ..., naming: ...}` の3つを持つ。
`grp` のような単一フィールドは使わない（v0 の設計ミス）。

---

## 5. 分子スキーマ（data/molecules/*.yaml）

`data/molecules/_TEMPLATE.yaml` を参照。要点のみ：

- `producers` / `receivers` は `{cell_id: {level, source, status, note}}`。
  level は `3`=主要 / `2`=明確 / `1`=文脈依存 / `0`=なし。
- `system_closure` は**書かない**。producers/receivers から build.py が自動判定する
  （`closed` / `out_only` / `in_only` / `external`）。手で書くと必ず食い違う。
- `own_data` は本人の実測のみ。日付・条件・値・単位を構造化して持つ。ここに一般知識を書かない。
- `effects` は `{fibrosis, inflammation, lipid, differentiation}` を -2〜+2。
  文献間で不一致なら `0` にして `note` に不一致の内容を書く。平均を取らない。

---

## 6. サイトの構造

```
/                  分子一覧＋5×5通信行列＋絞り込み
/molecule/IL6      分子ページ（受容体・経路・産生受容・測定法・自分の実測）
/pathway/il6-gp130 経路図（フェーズ2）
/theme/liver-inflammation  テーマ図（フェーズ3）
/families          3軸切り替えのファミリーツリー。「名前と実体のズレ」を一覧表示
/family/tgf-beta-superfamily  ファミリーページ（所属分子・分岐・例外）
/cell/ikc          細胞ページ（出す分子・受ける分子・分化元/先）
/todo              status: todo と inferred の一覧（作業リストとして使う）
```

`/todo` はフェーズ1から作る。**未検証がどこに残っているかが常に見える状態**にしておくことが、
このプロジェクトが信用できるものになるかどうかの分かれ目。

---

## 7. 作業の進め方（ユーザーとのやりとり）

ユーザーは要件定義・コーディングを専門としない。以下を守る。
- 「IL-15 を追加して」と言われたら、YAML を作り、`inferred` で埋め、出典が取れたものだけ `verified` にし、
  **何を verified にして何を inferred のまま残したかを報告する**。
- 大きな構造変更をする前に、影響範囲（何ファイルが変わるか）を先に伝える。
- 一度に4フェーズを並走させない。ユーザーが横道に逸れたら、フェーズ1が未完である旨を伝える。
