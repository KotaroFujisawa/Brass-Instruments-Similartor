# 参考文献とモデルでの用途

[English](../references.md) · 日本語

## 楽器の音響と形状

- Murray Campbell、Joël Gilbert、Arnold Myers（2021）、*The Science of Brass Instruments*、Springer。管体形状、マウスピースの共鳴、非線形伝搬、放射の近似の主な出典です。Yorkトランペットの管形状は図7.52（冊子373頁）から読み取りました。Conn 8Hの管体と5ALのカップは図4.22と図4.32bから再構成し、インピーダンスの比較には図4.46を用いました。図7.45と図7.46でAlexander C管チューバの実機と管形状を確認し、図6.9のフレアを使ってBB♭管の再構成に必要なベル形状を補いました。 [出版社・DOI](https://doi.org/10.1007/978-3-030-55686-0)
- Miranda Jackson、Gary Scavone（2025）、“A study of impedance of brass instruments and mouthpieces,” *Proceedings of Meetings on Acoustics* **58**、035019。Bach 5Cの形状データは図2の実測輪郭をデジタル化したものです。掲載図から輪郭座標を取り出し、単位をメートルに換算して、軸方向の原点を入口面に置きました。出典の論文はCC BY-NC 4.0です。 [DOI](https://doi.org/10.1121/2.0002165) · [ライセンス](https://creativecommons.org/licenses/by-nc/4.0/)
- Vennture Mouthpieces, Inc.（2019）、*VennCAD User Manual*、v2.0.3、16頁。資料の参照モデル例から、Yamaha BB-67C4の再構成に使う寸法と容積を取得しました。スロートの長さを仮定し、バックボア出口をチューバのマウスピース受け口に合わせた上で、これらの計測値に合う形状を求めています。 [マニュアル](https://www.vennture.mp/downloads/Documentation/VennCADHelp_v2-0-3.pdf)
- Antoine Chaigne、Jean Kergomard（2016）、*Acoustics of Musical Instruments*、Springer。7.3.4節など、マウスピースと管体を分布定数系として扱う音響モデルの背景です。 [出版社・DOI](https://doi.org/10.1007/978-1-4939-3679-3)
- Neville H. Fletcher、Thomas D. Rossing（1998）、*The Physics of Musical Instruments*、第2版、Springer。唇で駆動する金管楽器と、その励振の仕組みの背景です。 [出版社・DOI](https://doi.org/10.1007/978-0-387-21603-4)

## 唇の運動

- Seiji Adachi、Masa-aki Sato（1996）、“Trumpet sound simulation using a two-dimensional lip vibration model,” *The Journal of the Acoustical Society of America* **99**、1200–1209。唇先端の軸方向と横方向の運動、唇の接触、隙間の流体慣性、唇の移動に伴う体積流量のモデルの基礎です。このリポジトリの高音域フォルテ設定では、論文で調べた吹奏圧と唇の固有振動数の範囲を拡張しています。 [DOI](https://doi.org/10.1121/1.414601)
- Vincent Fréour、Christophe Vergez、Bruno Cochelin、Keita Arimoto（2025）、“Modelling the lips of an artificial trumpet player system,” *Proceedings of Meetings on Acoustics* **58**、035027。von Kármán板の定式化から、三次剛性を含む人工唇の低次元モデルを導いています。このリポジトリの三次項は学会要旨を参考にしています。係数の値とAdachi–Satoモデルとの組み合わせは実装で選んだものです。 [学会プログラム・要旨](https://ismra2025.org/wp-content/uploads/2025/05/ISMRA-2025-Final-1.3.pdf)
- Ukyo Inomata、Naoto Wakatsuki、Tadashi Ebihara、Yuka Maeda（2026）、“Numerical Investigation of Blowing-Pressure Effects on Artificial Lip Vibration in Brass Instruments,” *Journal of Signal Processing* **30**(4)、119–122。吹奏圧下での準静的な変形と幾何剛性が、人工唇の固有振動数に及ぼす影響を調べています。変位に依存する剛性を考える背景です。 [DOI・論文](https://doi.org/10.2299/jsp.30.119)

## 管内伝搬、損失、放射

- Harold Berjamin、Bruno Lombard、Christophe Vergez、Emmanuel Cottanceau（2017）、“Time-domain numerical modeling of brass instruments including nonlinear wave propagation, viscothermal losses, and lips vibration,” *Acta Acustica united with Acustica* **103**(1)、117–131。CPU版の双方向非線形伝搬、拡散表現による壁面損失、時間領域のフィードバックの参考です。 [著者原稿](https://arxiv.org/abs/1511.04247) · [DOI](https://doi.org/10.3813/AAA.919038)
- M. Jackson、Gary Scavone（2024）、“Modeled and measured impedance of brass instruments, mouthpieces, and bells,” *Proceedings of Meetings on Acoustics* **54**、035004。ベッセル関数で表す粘性・熱損失を使い、特にマウスピースの細いスロートに着目しています。マウスピース、管体、放射のモデルを実測と比較しており、局所的な円管の閉包を採用し、マウスピースの形状を詳しく扱う際の参考です。 [著者原稿](https://caml.music.mcgill.ca/lib/exe/fetch.php?media=publications%3Ajackson_poma_2024.pdf)
- Filipe Soares、Bruno Cochelin、Vincent Fréour、Christophe Vergez、Keita Arimoto（2025）、“Nonlinear acoustics in brass instruments using 2D complex modal bases,” *Acta Acustica* **9**、68。非線形伝搬、粘性・熱損失、外部放射音場を含む軸対称の低次元モデルを示しています。このプロジェクトでは、管内部に有限体積法を使い、等価球のモードから別途求めた放射負荷を結合しています。 [DOI](https://doi.org/10.1051/aacus/2025049)
- Leslie Greengard、Thomas Hagstrom、Shidong Jiang（2013年のプレプリント）、“The solution of the scalar wave equation in the exterior of a sphere.” 外向き球面波の解と時間領域での表現の参考です。 [著者原稿](https://arxiv.org/abs/1308.0643)
- Pauline Eveno、Jean-Pierre Dalmont、René Caussé、Joël Gilbert（2012）、“Wave Propagation and Radiation in a Horn: Comparisons Between Models and Measurements,” *Acta Acustica united with Acustica* **98**(1)、158–165。平面波と球面波の伝搬・放射の近似を、トランペットとトロンボーンのベルの実測と比較しています。 [DOI](https://doi.org/10.3813/AAA.918501)
- Fabrice Silva、Philippe Guillemain、Jean Kergomard、Bastien Mallaroni、Andrew N. Norris（2009）、“Approximation formulae for the acoustic radiation impedance of a cylindrical pipe,” *Journal of Sound and Vibration* **322**(1–2)、255–263。表IIにある、フランジのない開口の因果的な反射有理近似の出典です。使用する係数は0.167、1.393、0.457です。著者プレプリントでは表題が短くなっています。 [著者原稿](https://arxiv.org/abs/0811.3625) · [DOI](https://doi.org/10.1016/j.jsv.2008.11.008)

## 形状データの出典

`data/` のJSONには、初期のモデルで使った書籍のページ、式、デジタル化の手法、推定した寸法、仮定を記録しています。 [初期トランペットの再構成](../../data/sources.json) · [入射波を指定した非線形伝搬](../../data/nonlinear_sources.json) · [双方向CPUモデル](../../data/sources_bidirectional.json)

トロンボーンでは、書籍のグラフから読み取った形状を、メーカー公称のConn 8Hの管径0.547インチとDenis Wick 5ALの寸法で補いました。再構成したマウスピースは、書籍に示された全長、容積、共鳴を拘束条件として校正しています。 [Conn 8H](https://connselmer.com/instruments/trombones/conn-symphony-tenor-trombone-in-bb-8h) · [Denis Wickのトロンボーン用マウスピース](https://www.deniswick.com/trombone)

BB♭チューバのバルブ径20.5 mmとベル径450 mmには、AlexanderのModel 163の仕様を使いました。出典の管形状はエディンバラ大学楽器コレクションのAlexander C管チューバです。再構成では、管の中心線に沿った長さをBB♭管に換算し、共鳴間隔を校正しました。 [Alexander Model 163](https://gebr-alexander.de/products/b-oder-c-tuba-gebr-alexander-163m-messing-poliert-weite-bauart)
