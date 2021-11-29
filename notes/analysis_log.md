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

## 2021-08-19 TH
- Updated Draw2D_av.C: Added y histogram
- Made some corrections on Pol_power.C. Major editions are,
 - around l. 143: defined q_min: the minium value of the bins
 - around l. 183: defined twopirad2. The same angle as hq was used for hq2. Defined twopirad2 for hq2. Used around l. 232
 - around l. 230: division by the max kp value. In the previous code the same max kp as hq was used for hq2. Corrected tov kp2

## 2021-08-19  AH
- Made some corrections on Pol_power4.C. Major editions are,
 - around l. 143: defined q_min: the minium value of the bins
 - around l. 183: defined twopirad2. The same angle as hq was used for hq2. 
 Defined twopirad2 for hq2. Used around l. 232
 - around l. 230: division by the max kp value. In the previous code the same max kp as hq was used for hq2. Corrected tov kp2
 - Double_t nbin_q  = 300->60

 ## 2021-08-20  AH
 - Copy Pol_power4.C->Pol_power5.C . 

## 2021-08-20 TH
- Edited Pol_Power5.C. Fixed the problem that the y axis could not be set to log.
(hq[i]->GetYaxis()->SetRangeUser(1E-3,2.);) 
 - The section from l.299: cleared the canvas before plotting hq, hq2. Used if statement in order not to plot hq0 (Direct beam data)
 - From l.134: Adjusted the position of TLegend 

  ## 2021-08-20  AH
 - Added Pol_Power_30nm.C, Pol_Power_50nm.C, Pol_Power_90nm.C
 edit namestr, degstr, angle
 - Added q025_30nm.C (AFP magnetic field - Polarizing Power)

   ## 2021-08-21  AH
 - Edit q025_30nm.C (Added "" "")
 double nbin_qq=(qq-""q_min"")*nbin_q/(q_max-q_min);
 - Edit Pol_Power_30nm_temp.C (Added "" "")
    ibin_pol[i]= Int_t((q_cuts[i]-""q_min"")*nbin_q/(q_max-q_min));
    leg4->AddEntry(gr[i],Form("q=%.3f nm^{-1}",""q_min"" + (q_max-q_min)*ibin_pol[i]/""nbin_q""),"p");
 - Added  Pol_Power_50nm.C, Pol_Power_90nm.C
 Draw Graph mT vs Polarizing Power (lower right)
 - Added overview_14184125, overview_14193654, overview_17040703 
 (to compare 3 direct data)

 - Added AFP_RF_Pol.C, AFP_RF_Pol_deg.C
 (Draw AFP_RF mV vs Polarizing Power )
 - Added overview.C (to compare 3 direct data)

  ## 2021-08-31  AH

 - Added Draw2D_av_back.C to find background contribution
 - Added Pol_Power_30nm_scan_fit.C 
 - Added Pol_Power_50nm_scan_fit.C 
 - Added Pol_Power_90nm_scan_fit.C 

  ## 2021-09-08  TH
  - Exported individual root files from the scan 20210717002421. This is a scan without SF. This seems to have been stopped manually during cycle at index=9. The current applied at each cycle was checked from the log in 202107013_Polarizer/scan.txt (around L.4759). It actually looks it was stopped during I=0.47A (corresponding to index=10)->??. For now to be on conservative side, I exported data up to index=8. 
  - The corresponding magnetic field was calculated and placed as magnetic_20210717030252.csv.

  ## 2021-09-21  AH
  - Add P90nm_BG8020.C to write graph 80% up +20% down  
  - Add P90nm_BG10.C,P90nm_BG10.C to write B_RdataOFF_qq225_90nm_scan_fine_3.csv
  - Add P90nm_single.C to only write theoretical curve
  - Add all_B_90nm.C to write applied B vs R (Used for slides in JPS)
  - ~BG Background is taken into account

  ## 2021-10-05  AH
  - Add M1_q_R_I_trans_ref_4.C to write x-I,lambda-I,q-R,q-T graph throgh with M1 angular scan 

  ## 2021-10-24  AH
  - M1_pol_potential.C 
  model M1 polarization
  RPMT correction is not considered.!!(from Ichikawa)

  - fit_check_lambda.C
  lamda_vs_R のhistを描きたかったが失敗

  - M1_qc_m2_det_lam.C
  入射角度に依存しないモデル化　(4pi/lambda)
  
  - M1_lam_correct.C 
  After meeting, 変数をラムダにした関数でフィットするプログラム　
  パラメータ　λc(lmc) , m2(mm2)
  誤差も修正済み

  - M1_lam_correct1.C
  tanhxで減衰する横橋修論の式(f1)を用いたフィットと、q^-4で減衰する式(f0)でフィットした

  - fit_check_Spline_global2.C
  global fit to extract "saturation magnetization B(T)", V_Fe, sample thickness




