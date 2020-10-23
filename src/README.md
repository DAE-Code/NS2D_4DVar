# インプットファイルの説明（4DVar.inp）

インプットファイル4DVar.inpの中身は以下のようになっています．
```
2                  !--- Mode (1:TLM&ADJ Check, 2:4DV, 3:4DV restart)
1                  !--- Problem 1:Karman vortex, 2:Vortex advection
17 50              !--- Number of 4DV cycles, iteration on each cycle
0.0                !--- Error variance of pseudo measurement
0.1                !--- Measurement error variance
0.0                !--- Coefficient of background term
4 4 40             !--- Every iskip & jskip & tskip for measuremnt 
-2.6d0 6.0d0       !--- Left and right of measurement area
-2.6d0 2.6d0       !--- Bottom and top of measurement area
------------------------------------------------------------------------------
150 150            !--- DA window step, Total number of time step
100                !--- Reynolds number 
80                 !--- Number of mesh: jmax (imax x jmax, imax = 3*jmax)
0.02               !--- time step
1                  !--- 0:initial, 1:restart
1                  !--- 0:No output to screen,1:otherwise
10                 !--- Output plot3d interval
20                 !--- Output history interval
```

# ソースコードの説明

## NS2D_4DVar.f90

メインルーチンです．インプットファイル（4DVar.inp）の読み込みと変数の初期化サブルーチン呼び出し，そして，実行モードの選択（線形コード・アジョイントコードの確認か4次元変分法の実行）が行われます．

## mod_variables.f90

使用する変数を定義しています．

## sub_4dvar.f90

4次元変分法を実行するサブルーチンです．

## sub_check_foa.f90

線形コード・アジョイントコードを確認するためのサブルーチンです．

## m_random3.f90

平均0，分散1の正規乱数ベクトルを出力するサブルーチンです．変数の数n_varとアンサンブル数n_ensを入力して，そのサイズの正規乱数ベクトルを得ます．

## makefile

ソースコードをコンパイル・リンクして，実行ファイル4dvarを作成します．通常はsrcフォルダ内で`make`と入力します．デバッグモードでコンパイルする場合には`make debug`とします．

## sub_bc_outer.f90

計算領域の外側境界に流れの流入・流出などの境界条件を与えます．通常の2次元NSコードの外部境界条件に加えて，その線形コードおよびアジョイントコードが含まれています．

## sub_bc_wall.f90

無限長さ角柱（2次元計算領域における正方形）まわりの境界条件を埋め込み境界法によって与えるサブルーチンです．線形コードおよびアジョイントコードが含まれています．

## sub_hsmac_adj.f90

`sub_hsmac_fwd.f90`のアジョイントコードです．

## sub_hsmac_fwd.f90

圧力および速度の同時過緩和をHSMAC法によって行います．

## sub_hsmac_tlm.f90

`sub_hsmac_fwd.f90`の線形コードです．

## sub_initial.f90

使用する配列の確保，計算領域の定義を行います．`sub_addvtx(xc,yc)`は計算領域内にBurnham-Hallock渦を1つ設置するサブルーチンです．リスタートファイルの入力を行うサブルーチンも含まれています．

## sub_measure.f90

観測モデルを定義しています．どの範囲の計算格子点を何点おきに計測点とするかを定義してます．

## sub_plot3d.f90

流れ場の可視化用にPlot3dファイル形式で結果を出力します．Plot3dファイル（計算格子：mesh.g, 結果ファイル：*.q）はTecplotやParaviewで可視化できます．

## sub_rhs3rd_adj.f90

`sub_rhs3rd_fwd.f90`のアジョイントコードです．

## sub_rhs3rd_fwd.f90

流れ場の時間発展をオイラー陽解法（1次精度陽解法）で行います．対流項はKuwahara-Kawamuraスキーム，粘性項は2次精度中心差分で離散化します．

## sub_rhs3rd_tlm.f90

`sub_rhs3rd_fwd.f90`の線形コードです．
