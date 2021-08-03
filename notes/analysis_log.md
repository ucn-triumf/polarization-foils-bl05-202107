# analysis_log.md

## 2021-07-23 TH
Created the list of run # and scan logs in codes/

## 2021-07-27 AH
copy files "210713_FeSi  bin  Ichikawa" from "tools" to "codes"
change path "AFP.C  FeMirror.C  FeMirror2.C  Polarizer.C" 

## 2021-07-30 07:58 TH
made a simple code to draw a 2D histogram of the (x y) beam spot 

## 2021-07-31  AH
主に編集したのは、AFP2.C(色々変えてコード理解)
疑問点には?をつけた
その他疑問点
・q-Rの図でなぜ反射率>1となってしまうのかがわからない
答え　cut_refをcut_dirへ 150ぎょうめ

## 2021-08-2  AH
simple/FeMirror2.C　が正しいq-Rのスケーリング
FeMirror_lambda_counts.C　は一つのグラフのみプロット



|index|run #|Lav I(A)|real I(A) |mag B (kitaguchi)|AFP ON/OFF|refelction wave x |
20210715075452 	0	-0.0041	-0.32198	off	47.2
20210715081447 	0	-0.0041	-0.32198	off	47.09
20210715084835 	0.15	0.129145	-0.90759	off	47.04
20210715085349 	0.264	0.230411	-1.35266	off	47.2
20210715082606 	0.378	0.331677	-1.79772	off	47.2
20210715083711 	0.6	0.52888	-2.66443	off	47.2
20210715072653 	1.97	1.745851	-8.01302	off	47.18
20210715080233 	0	-0.0041	-0.32198	on	47.07
20210715082018 	0	-0.0041	-0.32198	on	47.13
20210715085144 	0.15	0.129145	-0.90759	on	47.1
20210715085714 	0.264	0.230411	-1.35266	on	47.19
20210715083141 	0.378	0.331677	-1.79772	on	47.24
20210715084052 	0.6	0.52888	-2.66443	on	47.21
20210715073913 	1.97	1.745851	-8.01302	on	47.11


## 2021-08-3  AH
正確な距離 20m?
angleの70.5, 56,51,64.71 がどこから来たか不明
Draw2D_av.Cの使い方
反射波のxの位置
