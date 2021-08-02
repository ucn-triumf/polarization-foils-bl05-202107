# analysis_log.md

## 2021-07-23 TH
Created the list of run # and scan logs in codes/

## 2021-07-27 AH
copy files "210713_FeSi, bin, Ichikawa" from "tools" to "codes"
change path "AFP.C, FeMirror.C, FeMirror2.C, Polarizer.C" 

## 2021-07-30 07:58 TH
made a simple code to draw a 2D histogram of the (x,y) beam spot 

## 2021-07-31  TH
主に編集したのは、AFP2.C(色々変えてコード理解)
疑問点には?をつけた
その他疑問点
・q-Rの図でなぜ反射率>1となってしまうのかがわからない
答え　cut_refをcut_dirへ 150ぎょうめ

simple/FeMirror2.C　が正しいq-Rのスケーリング

## 2021-08-02 TH
Added codes/simple/Draw2D_av.C
This code can be used to extract the peak value for a given root file 



