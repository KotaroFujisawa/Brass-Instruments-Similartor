# F/B♭ダブルホルンの軸対称モデル

実音B♭1からB♭5までのB♭メジャー上行音階を、フォルテでアタックを強くした吹奏条件で計算します。全29音で、最初の28音は各1.08秒、最後のB♭5は1.65秒です。各音の後ろに0.09秒の無音を入れるため、全体は34.50秒です。[公開WAV](../../samples/horn_Bb1_Bb5_major_f_strong_attack.wav) · [測定結果](../results/horn_verified_measurements.json)

## 管形状とマウスピース

B♭側の基準となる内径は、Alexander Mainz製ダブルホルン（1950年頃、エディンバラ大学所蔵EU1804）の実測を示す *The Science of Brass Instruments* の図7.34から読み取りました。図から読み取れる管径はベル途中の約66 mmまでです。そこから先は接続点で傾きの連続するベッセル形状（指数1.5）で補い、終端には現行Alexander 103の公称ベル径310 mmを使いました。[B♭側の管径データ](../../data/horn_bb_bore_fig734.csv) · [教科書](https://doi.org/10.1007/978-3-030-55686-0) · [楽器目録](https://www.euchmi.ed.ac.uk/ubl/ublh06.html) · [ベッセル形状によるベルの近似](https://era.ed.ac.uk/items/4800ccdb-f802-4910-8e30-8464851241c7) · [Alexander 103の仕様](https://gebr-alexander.de/en/products/f-b-horn-gebr-alexander-modell-103-vater)

F側は、共通のリードパイプとベルの間でB♭側の中央の円筒に近い部分を約0.985 m延長して構成します。追加長はF1とB♭1の半波長の差です。マウスピースを含む開放時の全長はB♭側が約2.88 m、F側が約3.87 mです。第1、第2、第3バルブはそれぞれ2、1、3半音分の管長を加える独立した管で表し、F側とB♭側でそれぞれの長さを設定します。演奏時は指定した側とバルブに対応する一続きの管を組み立て、曲管を中心線に沿って展開します。この軸対称化では、曲がり、バルブの分岐や押し替え動作、管壁の振動は省略しています。

マウスピースはAlexander 8Mの公称カップ内径17.5 mm、深さ28 mm、スロート内径4.0 mmを使います。カップの曲線は教科書の図7.55(c)の内側に湾曲する形状を参考にして、公称寸法に合わせました。スロート長5 mmとマウスピース全長66 mmを仮定し、バックボアは再構成しました。[マウスピースの管径データ](../../data/alexander_8m_mouthpiece_reconstructed.csv) · [Alexander 8Mの公称寸法](https://www.alexanderaustralia.com.au/product/alexander-french-horn-mouthpiece/) · [教科書](https://doi.org/10.1007/978-3-030-55686-0)

## 実音・記譜音と運指

F調ホルンの記譜音は実音より完全五度（7半音）高く、親指レバーで管を切り替えても記譜は同じです。実音B♭1は記譜F2、実音B♭5は記譜F6です。[演奏プリセット](../../data/performance_horn.json) の `note` と `written_note` にその対応を、`side`、`valves`、`partial` に使用する管、押すバルブ、選んだ倍音次数を記録しています。

通常の音域では運指表を参照し、低音では選択した倍音とバルブで決まる管長を照合しています。最低音のB♭1はB♭側の開放ペダル音です。高音では、例えば実音G5、A5、B♭5をB♭側の第16倍音の候補として扱い、標準運指表の範囲を上へ拡張しています。音ごとの管長補正、唇の固有振動数、吹奏圧も同じプリセットに記録しています。管長補正は、再構成した管体で目標の音程に合わせるためのモデル上の有効長補正です。[ヤマハの運指表](https://www.yamaha.com/en/musical_instrument_guide/horn/play/play003.html) · [International Horn Societyの低音運指](https://www.hornsociety.org/hornzone/301-fingerings-in-the-low-range) · [Alfredの運指表](https://images.alfred.com/pdf/BOP/FingeringCharts/French_Horn.pdf)

F側はC2〜F2とC3〜C4の12音、B♭側はB♭1、G2〜B♭2、D4〜B♭5の17音に使いました。運指はB♭メジャー音階上の音に対応しています。

## 管内の非線形伝搬と唇の連成

中心線に沿う座標を `z`、半径座標を `r` とし、管内を圧縮性Euler方程式で解きます。`rho` は密度、`u` と `v` は軸方向と半径方向の速度、`p` は圧力、`E` は単位体積あたりの全エネルギーです。

```text
U   = (rho, rho*u, rho*v, E)
F_z = (rho*u, rho*u*u+p, rho*u*v, (E+p)*u)
F_r = (rho*v, rho*v*u, rho*v*v+p, (E+p)*v)
p   = (gamma-1) * [E-rho*(u*u+v*v)/2]

d(r*U)/dt + d(r*F_z)/dz + d(r*F_r)/dr = (0, 0, p, 0)
```

同心リングの有限体積上で、圧縮波の急峻化や衝撃波を扱える保存形の流束を使います。面の状態はMUSCLとMC制限関数で再構成し、通常はHLLC流束、必要に応じてHLLE流束で更新します。時間積分はSSP-RK2です。粘性と熱による壁面損失は、各セルの因果フィルタによる運動量緩和で近似します。[軸対称Eulerソルバーの詳細](AXISYMMETRIC_MODEL.md)

唇はAdachi–Satoの2自由度モデルを用い、前後・上下の変位、閉鎖時の接触、隙間を通る空気と唇の運動が押しのける空気を扱います。カップ断面積を `S_cup`、唇の総体積流量を `U_lip` とすると、断面で平均した管からの戻り特性量 `h` とカップの音圧 `p_cup` の関係は `p_cup = h + rho0*c0*U_lip/S_cup` です。この音圧と隙間の圧力が唇を駆動し、その変位と流量が次の管内音場を変えます。吹奏圧を6.5 msで立ち上げ、持続時の設定に対して一時的に44%を上乗せします。上乗せ分は28 msかけて余弦曲線に沿ってゼロへ戻ります。[Adachi and Sato, *Trumpet sound simulation using a two-dimensional lip vibration model*](https://doi.org/10.1121/1.414601)

## 右手を入れたベルの放射と反射

ホルンのベルに右手を入れると、ベルでの反射が変わります。通常の手の位置と深く塞いだ状態の測定例は、教科書の第4.5.6節と図4.75を参照しました。YoshikawaとNobaraの図4(b)にある通常の手の位置の近似形状では、手を置く位置で管の内径が84 mm、開口の径が50 mmです。この局所的な面積比 `alpha = (50/84)^2 ≈ 0.354` をベル終端の等価開口に転用します。外部領域は軸対称の等価球とし、手の形状や左右非対称の放射は等価開口の近似に含めます。[教科書](https://doi.org/10.1007/978-3-030-55686-0) · [Yoshikawa and Nobara, *Acoustical Modeling of Mutes for Brass Instruments*](https://doi.org/10.1007/978-3-319-47292-8_5)

ベル断面積を `S_b` とすると、等価開口の面積は `S_e = alpha*S_b` です。二つの面で音圧と体積流量を連続にするため、等価開口のモード速度は `u_e = u_b/alpha` になります。この接続では音響仕事率が一致します。`R_e` を等価開口の反射行列、`I` を単位行列とすると、管側の反射行列 `R_b` は次式で求まります。

```text
beta = (1-alpha)/(1+alpha)
R_b  = (I + beta*R_e)^(-1) (R_e + beta*I)
```

断面の速度を4モードに射影し、反射を66状態の因果系で計算します。反射はベル面の境界条件から管内へ戻り、非線形伝搬を経て唇を駆動します。同じ等価球解から、ベル軸上1 mの放射音圧を求めます。`alpha=1` では開放ベルの放射モデルに一致します。[等価球による外部放射モデル](../trumpet_bell_exterior.md)

半径151.5 mmの試験格子では、10 Hzから22 kHzまでの評価点で周波数領域の反射行列と因果系の応答の差は最大約 `3.2e-14`、反射行列の最大特異値は丸め誤差の範囲で1以下でした。この検証は境界モデルの数値的な整合性と受動性を確認するものです。Dellの実測では、手の置き方による入力インピーダンスの変化は約500 Hz以上で顕著でした。実楽器への適用精度は、こうした測定との比較で評価できます。[実際の手をかたどった模型による測定](https://research-repository.uwa.edu.au/en/publications/investigation-into-factors-affecting-the-acoustics-and-playabilit/)

## 計算結果

出力したWAVを再測定したところ、全29音の持続部で目標音程との差は最大13.9セント（C5）でした。各音の後ろに90 msの無音があり、クリッピングは0サンプルでした。ベル軸上1 mで計算した持続部の音圧レベルは、29音の中央値で91.9 dB SPLです。

## 再計算

NVIDIA GPUと対応ドライバー、CUDAランタイムコンパイラーのNVRTCを利用できる環境で、リポジトリのルートから実行します。[CUDAのセットアップ](../../README.md)

```sh
python -m pip install -r requirements.txt
python src/render_horn_scale.py
```

全音の計算と検証が通ると、`output/horn_axisymmetric/horn_Bb1_Bb5_major_f_strong_attack.wav`（48 kHz、16 bit、モノラル）と `output/horn_axisymmetric/measurements.json` が出力されます。測定値には各音の周波数、目標からの音程差、カップ内と放射音の音圧が含まれます。WAVでは全音に同じゲインを使います。

一部の音の確認には `--only` と `--probe` を使えます。次の例は最低音と最高音をそれぞれ0.55秒計算して、確認用のJSONを出力します。

```sh
python src/render_horn_scale.py --only Bb1 Bb5 --probe 0.55
```
