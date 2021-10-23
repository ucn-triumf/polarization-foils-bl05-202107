MakeRPMTCalibrate.Cは、
対象となるtarget_list.rootのTTree Tに較正済み変数を追加したTTree Tcalib
を作成してtarget_calib.rootに保存する。

20210319211841_analyze.rootには較正用の関数が保存されている。

追加される変数はすべてDouble型で
xcalib: x座標 (mm)
ycalib: y座標 (mm)
wabs:   TOFとRPMTの検出効率から導いたウェイト
wrel:   RPMT上での相対的な検出効率から導いたウェイト
wcalib: wabsとwrelの積
tcalib: TOFの折り返しを行ったあとのTOF (ms)
である。
TOFの折り返しは、たとえば10 (ms)で折り返したとすると、
生データで10 (ms)より小さい5 (ms)だったものはtcalibは45 (ms)となり、
10 (ms)より大きい15 (ms)だったものはそのままtcalibは15 (ms)となる。

使用法は以下の通り。
target_list.rootはTTree Tが保存されたファイルである。
この場合、_list.root部分を_calib.rootに置き換えたtarget_calib.rootファイルが新しく作成される。
MLF sourceからRPMTまでの位置が18.0 (m)、
TOFの折り返しを10.0 (ms)で行うとき、
$ root -l -b -q 'MakeRPMTCalibrate.C+("20210319211841_analyze.root","target_list.root",18.0,10.0)'
で実行する。
MLF sourceからRPMTまでの位置とTOF折り返しは、wabsに影響する。
wabsを使わないとき（相対的なウェイトのwrelだけで良いとき）は適当な値を入れておけば良い。

target_calib.rootのTTree Tcalibを用いてDrawするときは、
TCalib->Draw("ycalib:xcalib","wcalib","colz");
とすれば、較正されたxy座標で検出効率ずみの2次元分布が得られる。
（座標がmm単位なのに注意）
検出効率の補正で、検出器上での相対的な効率のみで良い場合は
wcalibではなくwrelを使えば良い。
