# トランペットの物理シミュレーション

実音B♭3からB♭5まで、B♭メジャーの上行15音をフォルテの吹奏圧で計算する。各音の計算時間を0.72秒から1.08秒、最後のB♭5を1.10秒から1.65秒、音間と末尾の無音を0.06秒から0.09秒に延ばした。吹奏圧の立ち上がりとリリースの時間は固定した。完成WAVは18.12秒、48 kHz、16 bit、モノラル。各音は静止状態から個別に発振させ、全体に共通のゲインを適用する。

## 唇とマウスピース

Adachi–Sato (1996) の前後・上下変位、接触、隙間の非線形流量と慣性、唇運動が押しのける体積流量を使う。上下変位の静止位置からの差を `eta` として、復元力を `-k_y * eta * (1 + alpha * eta**2)` とした。係数 `alpha = 8e4 m^-2` は今回の設定値で、変位1 mmで三次項は線形項の8%になる。Fréourらが人工唇を対象にvon Kármán板のモデルから導いた三次剛性を、既存の二自由度モデルに加えた。二自由度系の係数は人体で同定していない。Inomataらの解析では、吹奏圧による準静的変形に伴う幾何剛性の増大が、人工唇の固有振動数を上げうると報告されている。

管の形状は *The Science of Brass Instruments* 図7.52のYork製B♭管トランペットから読み取った。近似マウスピースではカップ、直径3.7 mmのスロート、バックボアを別々に格子化した。バルブは等価的な管長の追加で表し、アンブシュアと有効管長を音ごとに調整した。吹奏圧は自励振動が得られるよう音域ごとに設定した。

## 管内の減衰とベルの境界

管内の非線形伝搬と反射は、軸対称圧縮性Euler方程式をMUSCL–HLLCの有限体積法で解く。急な波面を衝撃波として扱える保存形を使う。唇への戻り波は管内のその時刻の入口状態から `p' - rho0*c0*u_z` を面積平均して得る。カップ圧と唇流量をこの戻り波込みで求め、同じカップ圧を次の唇の更新に使う。

JacksonとScavoneの円管モデルと同じBessel関数の応答を、局所半径ごとに時間領域で評価する。以下で `s` はラプラス変数、`a` は管半径、`I_n` は変形Bessel関数、`nu` と `D` は空気の動粘度と熱拡散率である。

```text
w_v = s*a^2/nu
w_T = s*a^2/D
F(w) = 2*I_1(sqrt(w))/(sqrt(w)*I_0(sqrt(w)))
D_v = s*F(w_v)/(1-F(w_v))
D_T = (gamma-1)*s*F(w_T)
```

`D_v` は運動量、`D_T` は音圧の緩和レート。Bessel関数の漸化式で求めた応答を、減衰率と重みが非負の10項の緩和過程で各々近似し、粘性の直流成分 `8*nu/a^2` を加える。係数は局所的な円管の線形音響から求め、断面に分布させる。壁の熱・粘性境界層を空間格子で直接解像するには、異なる格子と方程式が必要になる。

ベルでは断面平均の進行圧力波に対してSilvaらの円管開口（フランジなし）の因果的な反射係数を用いる。`x=s*a/c0` とすると `R=-(1+0.167*x)/(1+1.393*x+0.457*x^2)`。他の半径方向成分は局所特性境界で流出させる。入口とベルで外から空気が流れ込む場合には基準状態のエントロピーを与え、流出時は管内状態のエントロピーを使う。小マッハ数の境界特性と局所の円管近似を使うため、ベル外の音圧は全体積流量の時間微分から求めた距離1 mの単極子近似である。Soaresらの軸対称外部音場まで含めた研究では、放射の指向性も直接扱っている。

## 計算と検証

NVIDIA GPU、対応ドライバとNVRTC、`requirements.txt` のPython依存ライブラリを用いる。Windowsで確認したNVRTCは `nvidia-cuda-nvrtc-cu12==12.4.127`。環境変数 `CUDA_NVRTC_DIR` でNVRTCのDLLディレクトリを指定できる。

```sh
python src/render_trumpet_modern.py
python src/verify_trumpet_modern.py
```

`measurements.json` には、目標音高、持続区間前半と後半の音程、音圧、周期性、境界の圧力残差を各音について記録した。`verification.json` は弱いパルスに対する開放運指の入力インピーダンス、およびB♭4の戻り波を有効・無効にした場合の比較を含む。インピーダンスの前後比較は同じ形状・格子の数値実験である。実機と照合するには、同じ個体・マウスピースの測定データが必要になる。

持続後半の音程の最大偏差は8.2セント（G4）。D4、E♭4、F4では波形に周期ごとの変動がある。弱いパルスで得た開放運指の最初の6共鳴は改修前後で0.7 Hz以内に収まり、第1ピークは113.8から96.8 MPa·s/m³へ低下した。B♭4では連成時のカップ音圧が5.37 kPa RMS、戻り波を切った比較計算では自励振動が消えた。

## 文献

- Adachi & Sato (1996), *Trumpet sound simulation using a two-dimensional lip vibration model*. [JASA 99, 1200–1209](https://doi.org/10.1121/1.414601)
- Fréour, Vergez, Cochelin & Arimoto (2025), *Modelling the lips of an artificial trumpet player system*. [ISMA-RR 2025 講演要旨](https://ismra2025.org/wp-content/uploads/2025/05/ISMRA-2025-Final-1.3.pdf)
- Inomata, Wakatsuki, Ebihara & Maeda (2026), *Numerical Investigation of Blowing-Pressure Effects on Artificial Lip Vibration in Brass Instruments*. [Journal of Signal Processing 30, 119–122](https://doi.org/10.2299/jsp.30.119)
- Jackson & Scavone (2024), *Modeled and measured impedance of brass instruments, mouthpieces, and bells*. [POMA 54, 035004](https://caml.music.mcgill.ca/lib/exe/fetch.php?media=publications%3Ajackson_poma_2024.pdf)
- Soares, Cochelin, Fréour, Vergez & Arimoto (2025), *Nonlinear acoustics in brass instruments using 2D complex modal bases*. [Acta Acustica 9, 68](https://doi.org/10.1051/aacus/2025049)
- Silva, Guillemain, Kergomard, Mallaroni & Norris (2009), *Approximation of the acoustic radiation impedance of a cylindrical pipe*. [著者稿・表II](https://arxiv.org/abs/0811.3625)
- Campbell, Gilbert & Myers (2021), *The Science of Brass Instruments*, 図7.52および第4章。
