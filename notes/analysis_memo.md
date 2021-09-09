# analysis_memo




## BL05 data structure


## summary of important measurements

### Direct measurement 
| run #           | SF   |I_LV(A)|B_kita(mT)| comments |
| --------------- | ---- |------ |--------- |--------- |
| 20210714184125  | OFF  | 1.97 | -8.01302 | probably  I=1.97A, SF:OFF |
| 20210714193654  | OFF  | 1.97 |-8.01302  | probably  I=1.97A, SF:OFF |
| 20210717040703  | OFF  |	2.0	 |-8.13014 | recorded after the sample c was taken out |

### Characterization of AFP-SF/measurement of beam polarization

These measurements were done with two mirrors (m1 and m2)

#### The common settings of the AFP-SF
- The DC current/voltage of the AFP-SF should have been kept unchanged. That is 3.270 V / 1.329 A (from photograph taken after the last day). The note on p. 50 of the logbook says  3.259 V / 1.33 A. So there might have been small fluctuations, but should be stable at 1 mA order.
- I could not find the frequency set on the RF coil of the AFP-SF. The resonance frequnecy of the circuit was 103.6 kHz. So the value close to that should hav been used (TH on 2021-08-21)

#### Relation of current and voltage of the RF coil

Characterization measruement recorded on p.37 of the logbook

|  V_set (mV)  |    I_monitor (A) |
| ---- | ---- |
|    100  | 0.3605 |
| 150 | 0.4913 |
| 500 | 1.380 |
| 800 | 2.13832 |
| 600 | 1.63504 |
| 760 | 2.03701 |
| 1000 | 2.64813 |



List of measurements with AFP-SF OFF (the last part of beam alignment)

| run # |theta_m2 (deg) |x_m2(mm) |comments |
| ----- | -------- | --------|------- |
| 20210714181508 |276.4 |78.0 | not explicitly written, but probably at I_LV=1.97A, B=-8.01302 mT |
| 20210714185238 |275.7 |78.0 | not explicitly written, but probably at I_LV=1.97A, B=-8.01302 mT |
| 20210714191741 |276.2 |78.0 | not explicitly written, but probably at I_LV=1.97A, B=-8.01302 mT |
|  | | |  |
|  | | |  |

List of measurements with AFP-SF ON

| run # |theta_m2 (deg) |x_m2(mm) | SF-RF (mV)|I_LV(A)	|B_kita(mT)	|comments |
| ----- | -------- | --------|------- |------- |--------|------- |
| 20210714204714 |275.7 |78.0 | 100 | 1.97 | -8.01302 | I=1.97A not explicitly written, but probably |
| 20210714205602 |275.7 |78.0 | 760 | 1.97 | -8.01302 | I=1.97A not explicitly written, but probably |
| 20210714210221 |275.7 |78.0 | 1000 | 1.97 | -8.01302 | I=1.97A not explicitly written, but probably |
| 20210714211037 |275.7 |78.0 | 500 | 1.97 | -8.01302 | I=1.97A not explicitly written, but probably |
| 20210714211642 |275.7 |78.0 | 300 | 1.97 | -8.01302 | I=1.97A not explicitly written, but probably |
| 20210714214337 |275.7 |78.0 | 760 | 0 | -0.32198 |  |
| 20210714215803 |276.2 |78.0 | 760 | 1.97 | -8.01302 | theta_m2 slightly changed |

### Manual scans of the sample Fe 30 nm

| run #          | I_labview (A) | real I(A) | mag B (kitaguchi) | SF ON/OFF | Reflection beam x (mm) |
| -------------- | -------------- | --------- | ----------------- | --------- | ---------------------- |
| 20210715075452 | 0              | -0.0041   | -0.32198          | off       | 47.2                   |
| 20210715081447 | 0              | -0.0041   | -0.32198          | off       | 47.09                  |
| 20210715084835 | 0.15           | 0.129145  | -0.90759          | off       | 47.04                  |
| 20210715085349 | 0.264          | 0.230411  | -1.35266          | off       | 47.2                   |
| 20210715082606 | 0.378          | 0.331677  | -1.79772          | off       | 47.2                   |
| 20210715083711 | 0.6            | 0.52888   | -2.66443          | off       | 47.2                   |
| 20210715072653 | 1.97           | 1.745851  | -8.01302          | off       | 47.18                  |
| 20210715080233 | 0              | -0.0041   | -0.32198          | on        | 47.07                  |
| 20210715082018 | 0              | -0.0041   | -0.32198          | on        | 47.13                  |
| 20210715085144 | 0.15           | 0.129145  | -0.90759          | on        | 47.1                   |
| 20210715085714 | 0.264          | 0.230411  | -1.35266          | on        | 47.19                  |
| 20210715083141 | 0.378          | 0.331677  | -1.79772          | on        | 47.24                  |
| 20210715084052 | 0.6            | 0.52888   | -2.66443          | on        | 47.21                  |
| 20210715073913 | 1.97           | 1.745851  | -8.01302          | on        | 47.11                  |
|                |                |           |                   |           |                        |

## Manual scans of the sample Fe 90 nm

| run #          | L_labview I(A) | real I(A) | H (mT)   | SF ON/OFF | Reflection beam x (mm) |
| -------------- | -------------- | --------- | -------- | --------- | ---------------------- |
| 20210716230107 | 0              | -0.0041   | -0.32198 | off       |                        |
| 20210716232122 | 2              | 1.7725    | -8.13014 | off       |                        |
| 20210716233530 | 2              | 1.7725    | -8.13014 | on        |                        |
| 20210717001854 | 0.37           | 0.324571  | -1.76649 | off?      |                        |
| 20210717004515 | 0              | -0.0041   | -0.32198 | off?      |                        |
| 20210717005023 | 0              | -0.0041   | -0.32198 | on        |                        |
| 20210717022140 | 0.265          | 0.2313    | -1.35656 | off       |                        |
| 20210717022920 | 0.265          | 0.2313    | -1.35656 | on        |                        |
|                |                |           |          |           |                        |

## Manual scans of the sample Fe 50 nm

| run #          | L_labview I(A) | real I(A) | H (mT)   | SF ON/OFF | Reflection beam x (mm) |
| -------------- | -------------- | --------- | -------- | --------- | ---------------------- |
| 20210717051355 | 0              | -0.0041   | -0.32198 | ON        |                        |
| 20210717052755 | 0              | -0.0041   | -0.32198 | OFF       |                        |
| 20210717054725 | 2              | 1.7725    | -8.13014 | OFF       |                        |
| 20210717055945 | 2              | 1.7725    | -8.13014 | ON        |                        |
|                |                |           |          |           |                        |





### List of automatic scans

| Scan name                                 | run #          | sample t | I_start (LV, A) | I_end (LV, A) | dI (LV, A) | SF ON/OFF | Comments                                                     |
| ----------------------------------------- | -------------- | -------- | --------------- | ------------- | ---------- | --------- | ------------------------------------------------------------ |
|                                           | 20210716193147 | 30 nm    | 0.15            | 0.26          | 0.01       | OFF       | Manually stopped due to beam trouble?                        |
| --                                        | 20210716201056 | 30 nm    | 0.15            | 0.26          | 0.01       | OFF       | Cannot identify the scan file, lost by overwriting?          |
| scan20210713_flipper_agilent_scan_1       | 20210716210153 | 30 nm    | 0.2             | 0.22          | 0.005      | ON/OFF    | Observed the transition between 0.15 and 0.22A               |
| scan20210713_flipper_agilent_scan_long_1  | 20210716220736 | 30 nm    | 0               | 2.0           | 2.0        | ON/OFF    |                                                              |
| scan20210713_flipper_agilent_scan_rough_1 | 20210716235619 | 90 nm    | 0.37            | 0.27          | 0.01       | OFF       | Wait: 80 s / kp 2000                                         |
| scan20210713_flipper_agilent_scan_rough_2 | 20210717002421 | 90 nm    | 0.37            | 0.47          | 0.01       | OFF       | Wait: 80 s / kp 2000, manually stopped during index=9        |
| scan20210713_flipper_agilent_scan_rough_3 | 20210717011121 | 90 nm    | 0.35            | 0.10          | 0.05       | OFF?      | Wait: 40 s / kp 1000                                         |
| scan20210713_flipper_agilent_scan_fine_1  | 20210717012326 | 90 nm    | 0.240           | 0.265         | 0.005      | ON/OFF    | The N was 13 (should be an even #), manually stopped? The scan log file lost |
| scan20210713_flipper_agilent_scan_rough_4 | 20210717024407 | 90 nm    | 0.265           | 0.325         | 0.01       | OFF?      |                                                              |
| scan20210713_flipper_agilent_scan_fine_3  | 20210717030252 | 90 nm    | 0.295           | 0.325         | 0.005      | ON/OFF    | Wait: 240 s / kp 6000                                        |
| scan20210713_flipper_agilent_scan_rough_5 | 20210717050342 | 50 nm    | 0.19            | 0.36          | 0.01       | OFF       | Wait: 20 s / kp 500                                          |
| scan20210713_flipper_agilent_scan_rough_6 | 20210717054010 | 50 nm    | 0.23            | 0.262         | 0.002      | OFF       | Creation time of the scan file is 05:21. Due to the beam stop? Wait: 20 s/kp 500 |
| scan20210713_flipper_agilent_scan_fine_5  | 20210717061701 | 50 nm    | 0.23            | 0.25          | 0.002      | ON/OFF    | Manually stopped at 07:25                                    |



### Scan export

- Use this command to export individual root file from scan root file. The resultant root files are saved under data_scans/ 

  >  sh  codes/cut_scans/run_makeelist_TH.sh 

- The scan log fie used for exporting is listed in the table below (also written in the bash script)

  | Fe thickness | Scan log                                     | Run#           | Comment                                                      |
  | ------------ | -------------------------------------------- | -------------- | ------------------------------------------------------------ |
  | 30 nm        | scan20210713_flipper_agilent_scan_1          | 20210716210153 |                                                              |
  | 30 nm        | scan20210713_flipper_agilent_scan_long_1     | 20210716220736 |                                                              |
  | 50 nm        | scan20210713_flipper_agilent_scan_fine_5_mod | 20210717061701 | The first end was missing -> complemented. The last start was removed |
  | 90 nm        | scan20210713_flipper_agilent_scan_fine_3     | 20210717030252 |                                                              |
  |              |                                              |                |                                                              |
  |              |                                              |                |                                                              |
  |              |                                              |                |                                                              |
  |              |                                              |                |                                                              |
  |              |                                              |                |                                                              |
  |              |                                              |                |                                                              |
  |              |                                              |                |                                                              |

  



