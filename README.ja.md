# 金管楽器シミュレーター

[English](README.md) · [日本語ドキュメント一覧](docs/ja/README.md) · [モデルの方程式と仮定](docs/ja/MODELS.md) · [参考文献](docs/ja/REFERENCES.md) · [形状データと演奏設定](data/README.ja.md)

B♭トランペット、テナー・トロンボーン、BB♭コントラバス・チューバの音を物理モデルで合成します。Adachi–Satoの2自由度唇モデルをマウスピースと管内音場に連成します。CUDAソルバーは、曲がった管を中心線に沿って展開した軸対称の管内で非線形な波の伝搬を解きます。ベルで放射された音を求めると同時に、管内へ戻る反射波が唇の振動に及ぼす影響も計算します。

## 音声サンプル

各サンプルは実音のB♭メジャー上行音階をフォルテで演奏します。基準音はA4 = 440 Hzです。22音を個別に静止状態から計算してつなぎ、各音を1.08秒、最後のB♭を1.65秒としました。各音の後に90 msの無音を入れ、合計26.31秒のWAVにしています。形式は48 kHz、16 bit、モノラルです。

| 音声 | 実音域 | 測定結果 |
| --- | --- | --- |
| [Yorkトランペット、Bach 5C](samples/trumpet_Bach5C_Bb3_Bb6_major_f_strong_attack.wav) | B♭3〜B♭6 | [音程と放射音圧](docs/results/trumpet_Bach5C_Bb3_Bb6_measurements.json) |
| [Conn 8Hトロンボーン、Denis Wick 5AL](samples/trombone_5AL_Bb2_Bb5_major_f_strong_attack.wav) | B♭2〜B♭5 | [基準のフォルテ](docs/results/trombone_5AL_Bb2_Bb5_measurements.json) |
| [トロンボーン、吹奏圧を上げた版](samples/trombone_5AL_Bb2_Bb5_major_louder_f_strong_attack.wav) | B♭2〜B♭5 | [持続時の吹奏圧を15%増加](docs/results/trombone_5AL_Bb2_Bb5_louder_f_measurements.json) |
| [トロンボーン、さらに吹奏圧とアタックを強めた版](samples/trombone_5AL_Bb2_Bb5_major_stronger_f_strong_attack.wav) | B♭2〜B♭5 | [音圧と発音時刻の比較](docs/results/trombone_5AL_Bb2_Bb5_stronger_f_measurements.json) |
| [再構成したBB♭チューバ、Yamaha BB-67C4](samples/BBb_tuba_Bb1_Bb4_major_f_strong_attack.wav) | B♭1〜B♭4 | [音程・音圧・発音時刻](docs/results/BBb_tuba_Bb1_Bb4_measurements.json) |

トロンボーンのアタックをさらに強めた版では、15%増圧の版から持続時の吹奏圧をさらに25%上げ、吹奏圧の立ち上がりを12 msから6.5 msに短縮しました。ベル正面1 mの放射音圧レベルは、22音の中央値で持続部が+3.8 dB、音が立ち上がってから最初の50 msが+6.0 dBでした。

音階のWAVはファイルごとに一つの出力ゲインを適用しています。2つの演奏に同じゲインを適用した抜粋で、吹奏圧とアタックを変えたときの音圧差を確認できます。音階全体を試聴するためのピークを抑えた版も収録しています。 [同じゲインで聴き比べ](samples/trombone_5AL_stronger_forte_before_after.wav) · [音階全体の試聴版](samples/trombone_5AL_Bb2_Bb5_major_stronger_f_listening.wav) · [以前の音源・再生成コマンド・チェックサム](samples/README.ja.md)

## インストール

開発にはPython 3.12を使用しました。NumPy、SciPy、Numba、Matplotlibが必要です。リポジトリのルートからコマンドを実行します。

```sh
git clone https://github.com/KotaroFujisawa/Brass-Instruments-Similartor.git
cd Brass-Instruments-Similartor
python -m venv .venv
```

Windows PowerShellでは `.venv\Scripts\Activate.ps1`、Linuxでは `source .venv/bin/activate` を実行して仮想環境を有効にします。その後、依存パッケージをインストールします。

```sh
python -m pip install -r requirements.txt
```

CPUモデルはここまでの設定で実行できます。軸対称モデルにはNVIDIA GPUとドライバー、およびCUDAランタイムコンパイラーのNVRTCが必要です。ソルバーはPythonの `ctypes` を通じてCUDA Driver APIを呼び出します。

```sh
python -m pip install nvidia-cuda-nvrtc-cu12==12.4.127
```

公開音源の計算にはNVIDIA RTX 2000 Ada GenerationとNVRTC 12.4.127を使用しました。Windowsでは `CUDA_NVRTC_DIR`、インストール済みのPythonパッケージ、`CUDA_PATH/bin` などからNVRTCを探します。Linuxでは `libnvrtc.so.12` と関連ライブラリをシステムローダーまたは `LD_LIBRARY_PATH` で参照可能にしてください。本プロジェクトでの動作確認環境はWindowsです。

## 3オクターブの音階を再生成する

```sh
python src/render_trumpet_bach5c.py
python src/render_trombone_detailed_scale.py
python src/render_trombone_detailed_scale.py --louder
python src/render_trombone_detailed_scale.py --stronger
python src/render_tuba_scale.py
```

WAVと音ごとの測定結果は、上のコマンド順に `output/trumpet_bach5c/`、`output/trombone_detailed_scale/`、`output/trombone_louder_forte/`、`output/trombone_stronger_forte/`、`output/tuba_axisymmetric/` に保存されます。トロンボーンの3つの版をこの順で実行すると、一つ前の版に対する音圧差も計算します。演奏設定には各音の唇の固有振動数、吹奏圧、静止時の唇の開き、運指またはスライド位置、管長補正が含まれます。

まず1音を短く試す場合は、次のコマンドで選んだ音の診断結果を保存できます。音階全体のWAVは上のコマンドで生成します。

```sh
python src/render_trumpet_bach5c.py --only Bb4 --probe 0.45
python src/render_trombone_detailed_scale.py --only Bb3 --probe 0.45
python src/render_tuba_scale.py --only Bb1 --probe 0.45
```

## CPUモデルを実行する

非線形な唇の振動と線形管体を連成したモデルは、次のコマンドで保存済み設定の音階と診断図を生成できます。

```sh
python src/render_scale.py --render
python src/verify.py
```

一次元波動ガイドでは、管内の波が双方向に非線形に伝わり、戻り波が唇に作用します。

```sh
python src/render_bidirectional.py
python src/verify_bidirectional.py
python src/verify_bidirectional.py --full
python src/plot_bidirectional.py
```

入射波を与える非線形伝搬モデルは `python src/render_nonlinear.py` で計算し、`python src/verify_nonlinear.py --full` で検証します。このモデルでは唇による励振を線形共鳴器で求め、その波形を非線形伝搬の入力にします。双方向波動ガイドとCUDAソルバーでは、戻り波も時々刻々の音場から求めます。

## モデルと形状

YorkトランペットとConn 8Hの管径は、Campbell、Gilbert、Myersの *The Science of Brass Instruments* に掲載された図から読み取りました。Bach 5Cには公表された実測断面を用い、Denis Wick 5ALはカップ断面を読み取って長さ・容積・単体共鳴から校正しました。BB♭チューバはAlexanderのC管の管径曲線をBB♭管に合わせ、Alexander 163の公称バルブ径とベル径を反映しています。Yamaha BB-67C4は、公開測定値の寸法と容積を基に再構成しました。 [形状データの出典](docs/ja/REFERENCES.md)

CUDAソルバーは圧縮性Euler方程式に基づき、密度、軸方向と半径方向の運動量、全エネルギーを時間積分します。同心リングの有限体積にMUSCL再構成とHLLC流束を用い、必要に応じてHLLE流束へ切り替えます。時間積分は2段のSSP Runge–Kutta法です。保存形の離散化で衝撃波を扱い、管から戻る音響特性量を用いて唇の流量とカップの圧力を結びます。

3オクターブ版は、円管の周波数依存の粘性・熱損失を考慮します。ベル開口の4モードを等価球の外部音場に結合し、反射波とベル正面1 mの放射音圧を求めます。トランペットでは流れに垂直な向きの唇の復元力に、三次の剛性項も加えています。軸対称化と等価球面は幾何形状の近似です。バルブ通路の曲がり、管壁の振動、乱流ジェット、奏者の声道、室内の反響を扱うには、モデルの拡張が必要です。 [方程式と仮定](docs/ja/MODELS.md) · [楽器ごとの技術ノート](docs/ja/README.md)

## 数値検証

共通のEulerソルバーではSodの衝撃波管、保存性、円筒管の半径方向モードを検証します。ベル外部モデルでは安定性、受動性、相反性を調べ、球面波の反射を別の方法でも計算して比較します。

```sh
python src/verify_axisymmetric.py
python src/verify_bell_exterior.py
```

完成したWAVを別に解析すると、持続部の音程誤差はBach 5Cトランペットで最大約12セント、5ALトロンボーンの3版で最大12セント以内、BB♭チューバで最大約10セントでした。弱い周期変動が見られた音は測定結果に記録しています。数値は再構成した形状と調整済みの演奏設定による合成音の結果です。 [測定と数値検証の一覧](docs/results/README.ja.md)

## リポジトリの構成

| パス | 内容 |
| --- | --- |
| `src/` | シミュレーション、音声生成、検証、作図 |
| `data/` | 管径とマウスピースの形状、演奏設定、出典 |
| `samples/` | 公開WAV、再生成手順、チェックサム |
| `docs/ja/MODELS.md` | 方程式、前提条件、検証範囲 |
| `docs/ja/REFERENCES.md` | 参考文献と各モデルとの関係 |
| `docs/results/` | 測定結果と数値検証 |
| `docs/figures/` | 診断図 |
| `docs/ja/` | 日本語の技術ノート |
| `output/` | ローカルで生成する音声、時系列、診断結果 |
| `tmp/` | ローカルで生成するキャッシュと作業ファイル |

## 参考文献

形状や音響特性にはCampbell、Gilbert、Myers（2021）、唇の2自由度モデルにはAdachi、Sato（1996）を用いています。人工唇の剛性、円管の粘性・熱損失、ベルの放射には、それぞれの研究を参照しました。書誌情報、図番号、データの利用条件は参考文献に列挙しています。 [参考文献の一覧](docs/ja/REFERENCES.md)
