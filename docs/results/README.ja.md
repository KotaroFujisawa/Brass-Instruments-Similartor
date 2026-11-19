# 公開測定結果

[English](README.md) · [音声サンプル](../../samples/README.ja.md)

公開した音声の計算値と解析結果です。音声生成や検証のスクリプトを再実行すると、`output/` に新しいレポートが出力されます。公開レポート内の音声ファイルのパスはリポジトリのルートからの相対パスです。

## 3オクターブの音階とモデルの検証

| レポート | 内容 |
| --- | --- |
| [Bach 5C を使ったトランペット](trumpet_Bach5C_Bb3_Bb6_measurements.json) | 22音の計算結果と、生成した PCM ファイルを別途解析した結果 |
| [トランペットの高音域比較データ](trumpet_high_measurements.json) | Bach 5C 版との比較に使った以前の B♭3–B♭6 の結果 |
| [Wick 5AL を使ったトロンボーン](trombone_5AL_Bb2_Bb5_measurements.json) | B♭2–B♭5 の各音の計算結果と、生成した WAV の独立した解析 |
| [5AL トロンボーンの音量を上げたフォルテ](trombone_5AL_Bb2_Bb5_louder_f_measurements.json) | 詳細形状を使った5ALの標準音階を基準とした、音圧レベルと PCM レベルの比較 |
| [5AL トロンボーンのアタックを強めた音階](trombone_5AL_Bb2_Bb5_stronger_f_measurements.json) | 発音時刻をそろえた区間と固定時刻の区間におけるアタックの比較、持続部のレベル、ピークを制限した試聴音源の測定値 |
| [Yamaha BB-67C4 を使ったBB♭チューバ](BBb_tuba_Bb1_Bb4_measurements.json) | B♭1–B♭4 の計算結果と PCM の独立した解析。発音のタイミングと発音初期の低音の変調を含みます。 |
| [トランペットのテーパーとベルの比較](bell_comparison.json) | 以前のテーパー、細分化したテーパー、連成した管外放射モデルを比較した結果。[テーパーの図](../figures/trumpet_taper_refinement.png) |
| [ベル境界の検証](bell_boundary_verification.json) | CUDA 計算のベル境界の圧力と、CPU で周波数領域から再構成した圧力の比較 |
| [ベルの管外放射の検証](bell_exterior_verification.json) | 開口面でのモードへの射影、状態空間モデルの近似、受動性、受音点への伝達特性の確認 |
| [トランペットのマウスピース比較](trumpet_mouthpiece_results.json) | 以前のマウスピース形状と Bach 5C の実測輪郭を用いたモデルで、形状、インピーダンス、吹奏時の結果を比較。[インピーダンスの図](../figures/mouthpiece_impedance.png) |

## 以前の結果

| レポート | 内容 |
| --- | --- |
| [トランペットの2オクターブ音階の測定](trumpet_modern_measurements.json) | 各音の音高、発音初期と持続部の音高推定、圧力、周期性、ソルバーの設定 |
| [トランペットの2オクターブ音階の検証](trumpet_modern_verification.json) | 弱いパルスによる入力インピーダンス、ベルの反射係数の大きさ、B♭4のフィードバックの確認 |
| [Euler ソルバーの試験](axisymmetric_solver_tests.json) | 静止状態、保存性、Sod 衝撃波管の格子細分化、円筒管の半径方向モードの確認 |
| [初期の軸対称モデルによるトランペット](trumpet_axisymmetric_measurements.json) | 以前の軸対称モデルによる音階の音ごとの結果 |
| [トロンボーンの測定](trombone_measurements.json) | 標準音階の結果、プリセット、音ごとの診断値 |
| [トロンボーンの形状検証](trombone_geometry_verification.json) | マウスピースの容積と共鳴、スライド長、インピーダンスのピーク |
| [トロンボーンの連成検証](trombone_coupling_tests.json) | B♭3のフィードバックと、格子を細分化した場合の検証 |
| [トロンボーンのアタックを強めた2オクターブ音階](trombone_strong_attack_measurements.json) | アタックの設定、発音時刻の比較、音高とレベルの変化 |
| [線形管体モデルによるトランペット](trumpet_linear_measurements.json) | 初期の CPU 版の音階 |
| [双方向モデルによる CPU 版トランペット](trumpet_bidirectional_measurements.json) | 非線形モデルと線形モデルを比較した指標 |
| [所与の波形を伝搬させたトランペット](trumpet_one_way_measurements.json) | 非線形伝搬の指標と、条件をそろえた基準音源との比較 |

記録した計算時間は、ハードウェアや計算時に並行して動いていた処理に依存します。`common_gain_fs_per_pa` や `gain_fs_per_pa` などのフィールドには、それぞれの音声生成時に PCM へ変換するために適用したゲインを記録しています。
