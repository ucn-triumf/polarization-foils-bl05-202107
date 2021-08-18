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


## 2021-08-03  AH
質問まとめ
正確な距離 20m?
angleの70.5, 56,51,64.71 がどこから来たか不明
Draw2D_av.Cの使い方
反射波のxの位置

## 2021-08-03 19:39 TH
- Added README
- Fixed the order of the table of run # on notes/run_list.md. Re-exported it by codes/python/parse.python
- Updated the geometrical infromation of codes/simple/q_R_I_off.C and codes/simple/q_R_I_on.C
- AH reported this morning that a behavior of (R,q) plot produced by q_R_I_off.C for data at i=1, with a legend "AFP OFF B=8.01 mT", but I didn;t observe the same. Bug in a local file on AH's PC?

## 2021-08-05  AH
- (mistaken) Added Pol_Power.C (only draw Graph "q - Polarizing Power" )
- (mistaken) Added Pol_Power2.C (draw Graph "q - Polarizing Power" with other 3 graphs)

## 2021-08-18  AH
- Added Pol_Power4.C (draw Graph "q - Polarizing Power", Lower left)
(upper left -> off, upper right -> on)

Pol_Power3.C -> Pol_Power4.C Correction Point
if(useMRfirst) kp2[i] = tup2[i]->GetMaximum("mp");//"mp2"->"mp"
    else kp2[i] = tup2[i]->GetMaximum("kp");//"kp2"->"kp"
