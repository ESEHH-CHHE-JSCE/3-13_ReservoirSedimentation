# 3-13_ReservoirSedimentation
1. プログラム概要
* ダム貯水池堆砂をSMAC法とCIP-CLS2法を組み合わせた1次元河床変動モデルで求めるプログラム．
* 詳細は，2024年度版水理公式集例題集の例題3.13を参照のこと

2. コンパイル方法
* ifort プログラム名.f90

3. 実行方法
* プログラム名.exe

4. 入出力データ
* 入力データはプログラム中に記述．入力データファイルは不要
* 出力データは下記のとおり
  * TEST-TIME_1D_06_S.dat : 3時間毎の全格子点の流速，水深，水位，河床位，浮遊砂濃度の出力
  * TEST-TIME_1D_06_S_IN.dat : 1分毎の流入境界端の掃流砂量，浮遊砂巻き上げ量，沈降速度×河床面近傍基準面濃度，浮遊砂濃度，水深，浮遊砂濃度×水深，流速，浮遊砂フラックスの出力
  * TEST-TIME_1D_06_S_OUT.dat : 1分毎の流出境界端の掃流砂量，浮遊砂巻き上げ量，沈降速度×河床面近傍基準面濃度，浮遊砂濃度，水深，浮遊砂濃度×水深，流速，浮遊砂フラックスの出力
