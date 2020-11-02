# NS2D_4DVar

2次元ナビエ・ストークス方程式を用いた角柱まわりのカルマン渦列流れにおいて，4次元変分法(4D-Var)によるデータ同化を行うためのプログラムです．「データ同化流体科学－流動現象のデジタルツイン－」（共立出版）の内容とリンクしておりますので，書籍と合わせてご利用頂くと効果的です．本プログラムは第5章の計算に利用したものです．

# Linux環境での実行方法

環境が構築済みの場合には以下のようにソースコードの入手およびコンパイルを行います．
```
$ git clone https://github.com/DAE-Code/NS2D_4DVar
$ cd ./NS2D_4DVar/src
$ make
```
計算結果はグラフやコンター図で確認できます．計算結果の可視化については[こちら](https://github.com/DAE-Code/NS2D_DataAssimilation)をご参照ください．計算からグラフ化までを一括して処理できるスクリプトを用意しています．

# 環境構築

プログラムのコンパイルにはFortranコンパイラおよびmakeを使用します．Redhat系Linux（CentOSなど）の場合には，`sudo yum install gfortran`，Debian系の場合（Ubuntuなど）には`sudo apt install gfortran`のようにしてFortranコンパイラgfortranをインストールすることができます．makeが入っていない場合にも同様にインストールしてください．
