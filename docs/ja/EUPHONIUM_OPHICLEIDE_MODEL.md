# ユーフォニアムとB♭管オフィクレイド

ユーフォニアムは実音B♭1〜B♭5、B♭管オフィクレイドは実音B♭1〜B♭4を対象とします。いずれもフォルテでアタックを強めたB♭メジャーの上行音階です。基準音はA4 = 440 Hzです。ユーフォニアムの3オクターブ版（B♭2〜B♭5）も再生成できます。

## ユーフォニアムの形状とバルブ

管体は *The Science of Brass Instruments* 図2.21のユーフォニアムの曲線を177点で読み取ったものです。管体の中心線長は約2.713 m、入口径は約10.6 mm、ベル径は約279 mmです。図には個体を特定できる情報がありません。 [書籍](https://doi.org/10.1007/978-3-030-55686-0) · [管体データ](../../data/euphonium_bore_fig221.csv)

マウスピースはBoosey & Co.の資料個体13.15の寸法を使いました。全長73.1 mm、リム内径26.9 mm、カップ深さ36.0 mm、スロート径6.1 mmで、シャンク内径は約9.8 mmです。内面はDenis Wick 5AL用の曲線を計測値に合わせて再構成し、スロートの一定径区間とレシーバーの接続部をそれぞれ4 mmと仮定します。 [Royal Birmingham Conservatoireの収蔵記録](https://www.bcu.ac.uk/conservatoire/research/hic/the-collection/mouthpieces/13-15) · [マウスピースの再構成曲線](../../data/boosey_euphonium_mouthpiece_reconstructed.csv)

4バルブのコンペンセイティング方式を近似します。第1〜4バルブの公称降下量は順に2、1、3、5半音です。第4バルブを使わないときは各管路の付加長を加えます。第4バルブを含むときは、その補償管を含む総延長を半音数の合計から求めます。付加長と音別の管長補正はバルブ区間に入れます。 [B♭管4バルブ機の収蔵例](https://www.euchmi.ed.ac.uk/ubl/ublh10.html)

## オフィクレイドの形状と指孔

管体はEdinburgh University Collection of Historic Musical InstrumentsのGautrot製B♭管バス・オフィクレイド（EU3590）の管径資料から再構成しました。*The Science of Brass Instruments* 図7.34とCampbell、Chick、Myers（2014）の図2を使い、ベル末端の形状は同系統の個体の直径215 mmを参照して補いました。管体の軸方向寸法は計算時に調律のため0.985倍します。 [EU3590の収蔵記録](https://www.euchmi.ed.ac.uk/ubl/ublh07.html) · [管径の測定図](https://www.historicbrass.org/images/hbj/hbj-2014/HBSJ_2014_JL01_005_Myers.pdf) · [比較個体のベル寸法](https://instrumantiq.com/pages_html/galeries_AVC-1114-ophi-Gautrot.htm) · [再構成した管体](../../data/ophicleide_gautrot_eu3590.csv)

マウスピースはSegressが転載したCaussinusの寸法を参照しました。全長72 mm、カップ径27 mmと深さ34 mm、スロート径7 mm、出口径10 mmを使い、内面の曲率、長さ4 mmの一定径スロート、長さ3 mmの接続部を近似しました。 [Segressの論文、24頁](https://digital.library.unt.edu/ark:/67531/metadc935706/)

管にはキーで開閉する側孔を9個配置します。第1孔は通常開いており、対応するキーを押すと閉じます。第2〜9孔はキーを押すと開きます。孔の位置と初期運指は教科書4.4節とB♭管の運指表を参照しました。孔径と管壁から開口端までの高さは推定値です。D2では資料の比較条件に合わせて第5孔だけを開け、C3では試算で安定した、第1孔と第4孔を開ける代替運指を採用しました。音ごとの管長補正は再構成モデルの調律値で、管体と孔位置の軸方向寸法に反映します。ベルまでの管路をすべての運指で計算領域に含めるため、開孔より下流にある管とベルの影響も残ります。 [B♭管オフィクレイドの運指表](https://www.serpentwebsite.com/o_finger_big.pdf)

## 管内の連成と放射

巻いた管の中心線を展開し、周方向に対称な有限体積格子を作ります。管内は粘性・熱損失を含む軸対称の圧縮性Euler方程式を使い、時間積分には衝撃波を扱える保存形を用います。Adachi–Satoの2自由度唇モデルは、吹奏圧とカップ内の圧力で駆動されます。唇の開閉と変位は管への流入を変え、反射波はカップの圧力を通じて唇に戻ります。ベルでは開口速度の4つの軸対称モードを等価球の外部音場に接続します。 [唇モデル](https://doi.org/10.1121/1.414601) · [管内ソルバー](AXISYMMETRIC_MODEL.md) · [等価球面のベル放射](../trumpet_bell_exterior.md)

オフィクレイドの開孔では、円形開口を軸方向に分割した面積を、対応する格子の外周壁面に配分します。側孔ごとに面積平均した圧力の特性量から、開口に垂直な流速を求めます。管壁から開口端までの短い部分（チムニー）の慣性と、放射による周波数依存の損失を負荷に含めます。各格子面の流束は、開口と残りの壁面の流束をそれぞれの面積に応じて足し合わせます。

孔半径を $a$、断面積を $S_h=\pi a^2$、物理長と入口の補正を足した有効チムニー長を $L$ とします。流量に対する音響インピーダンス $Z_h(s)$ は次式です。

$$
Z_h(s)=\frac{\rho_0 c}{S_h}\left[\frac{sL}{c}+
\frac{0.61^2}{0.25}\frac{sT/c}{1+sT/c}\right],
\qquad T=\frac{0.25}{0.61}a.
$$

入口の補正には、側孔と本管の接続部に生じる内部の端補正と、曲面との接続による体積の補正を含めます。第2項は周波数依存の放射負荷です。低周波での放射端補正 $0.61a$ と、規格化した放射抵抗 $(ka)^2/4$（$k$ は波数）を再現します。放射負荷を因果的な緩和過程で解き、開孔の圧力と流量を管内と毎時刻連成します。 [Lefebvre & Scavone（2012）式8、10、12](https://doi.org/10.1121/1.3685481)

出力音圧はベル正面1 mを基準にします。オフィクレイドでは、各開孔の体積流量の時間変化から自由空間の単極子放射も計算します。孔から受音点までの距離を展開した管の座標で近似し、距離減衰とベル放射に対する時間差を反映して加えます。閉じた孔は剛体壁とし、パッドのたわみと閉じたチムニー内の空気の圧縮は省略しています。孔の周方向位置による指向性や、複数の孔からの放射がほかの孔の入口圧力に及ぼす影響も簡略化しています。

## 側孔境界の検証

半径16 mmの直円管に、入口から0.41 mと0.76 mの位置に2つの孔を設けて弱いパルスを送りました。開口面積の和は設定値と一致し、静止状態の圧力と流速の変化は数値上ゼロでした。30 Paを基準とするパルスに対し、各孔での流量の立ち上がりは約1.510 msと2.531 msでした。開孔時と全閉時の入口圧力を比較すると、最初の反射は約2.755 msで現れ、伝搬距離からの概算である2.841 msと0.086 msの差でした。開孔を持たない管の衝撃波管、保存性、円筒管の半径方向モードの検証も行いました。

```sh
python src/verify_toneholes.py
python src/verify_axisymmetric.py
```

## 音階の再生成

公開音源を再生成するには、CUDAの実行環境でリポジトリのルートから実行します。

```sh
python src/render_euphonium_ophicleide.py euphonium --extended
python src/render_euphonium_ophicleide.py ophicleide
```

各音を静止状態から独立に計算します。6.5 msで吹奏圧を立ち上げ、標準設定では持続圧の45%を一時的に上乗せして28 msで戻します。音別の唇の設定、吹奏圧、運指、管長補正を参照できます。[ユーフォニアムの設定](../../data/performance_euphonium.json) · [ユーフォニアムの追加低音の設定](../../data/performance_euphonium_pedal.json) · [オフィクレイドの設定](../../data/performance_ophicleide.json)

通常の音は1.08秒、最後のB♭は1.65秒で、各音の後に90 msの無音を入れます。WAVは48 kHz、16 bit、モノラルです。音声と `measurements.json` は、ユーフォニアムが `output/euphonium_extended/`、オフィクレイドが `output/ophicleide_axisymmetric/` に出力されます。測定ファイルには音程、音圧、周期性、格子と時間刻みを記録し、オフィクレイドでは開孔と孔からの放射音圧も記録します。

| 楽器 | 実音域 | 音数 | 長さ | 公開ファイル |
| --- | --- | ---: | ---: | --- |
| ユーフォニアム | B♭1〜B♭5 | 29 | 34.50秒 | [WAV](../../samples/euphonium_Bb1_Bb5_major_f_strong_attack.wav) · [測定結果](../results/euphonium_Bb1_Bb5_verified_measurements.json) |
| B♭管オフィクレイド | B♭1〜B♭4 | 22 | 26.31秒 | [WAV](../../samples/ophicleide_Bb1_Bb4_major_f_strong_attack.wav) · [測定結果](../results/ophicleide_verified_measurements.json) |

短い試算は次のコマンドで実行できます。

```sh
python src/render_euphonium_ophicleide.py euphonium --only Bb2 F3 --probe 0.45
python src/render_euphonium_ophicleide.py ophicleide --only Bb1 F2 --probe 0.45
```

### ユーフォニアムの3オクターブ版

```sh
python src/render_euphonium_ophicleide.py euphonium
```

B♭2からB♭5までの22音を計算し、`output/euphonium_axisymmetric/` に音声と測定結果を出力します。公開した4オクターブ版は、この設定に開放運指のペダル音B♭1と、バルブを使うC2〜A2の6音を加えたものです。追加の低音7音にも同じ管体・マウスピース・放射モデルを使います。[補償式ユーフォニアムの運指表](https://olemiss.edu/lowbrass/lowbrassmethods/charts/euphoniumbcfingering.pdf)

## 計算結果

両楽器の3オクターブ版のWAVを独立に解析し、全22音、全長26.31秒、クリッピング0を確認しました。持続部の音程誤差は、ユーフォニアムで最大13.3セント、オフィクレイドで最大12.5セントでした。オフィクレイドは22音中19音で側孔を開き、いずれも開孔からの放射を確認しました。C3、E♭4、G4には弱い周期倍化が見られました。

ユーフォニアムの4オクターブ版でも、全29音、全長34.50秒、クリッピング0を確認しました。持続部の音程誤差は最大13.3セントで、追加した開放B♭1の持続周波数は約58.24 Hzです。B♭2〜B♭5のWAVのサンプル値は、3オクターブ版の対応する区間と完全に一致しました。C2とD2の発音直後には弱い周期倍化が見られました。
