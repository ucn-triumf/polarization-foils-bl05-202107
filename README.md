# polarization-foils-bl05-202107
 updated: 2021-08-03 Takashi Higuchi

## Description of the repository 
This is a repository for managing codes to analyze data of cold-neutron reflectometry measurement which took place during beamtime in July 2021 at J-PARC MLF BL05.

## Description of the repository structure 
The contents of the respoistory are organized as follows:

- codes/ : place for codes written by H. Akatsuka and T. Higuchi   
 
- data/ : the directory which contains raw data. Omitted for github sync by .gitignore
	- data/210713_SiFe : directory copied from niki, the linux machine that acquired data from the NIM module. Contains raw RPMT data.
	- data/202107013_Polarizer : directory copied from kirk, the winodws computer which was used for controlling devices 
	- data/gatenetlog : directory copied from niki. Contains scan logs.
	- data/sample_magenticfield_0715 : results of magnetic field measurements by a Hall probe on the July 15th after the first two days of the beamtime. It was done while a sample was in place in the magnet  
	- data/sample_magenticfield_0717 : results of magnetic field measurements by a Hall probe on the July 17th after the beamtime. The sample was removed and the probe was placed near the center of the magnet.

- docs/ : place for documentations
	- docs/AFP_SF.pdf : summary of characterization measurements of AFP spin flipper  	
	- docs/BL05v21_17_82.pdf : scan of the BL05 log book  	
	- docs/electromagnet_drawing.pdf : the drawing of the electromagnet used for the measurement.
	- docs/I_correction1.pdf : H. Akatsuka's summary on calibration of currents applied to the magnet and the magnetic field produced
	- docs/Measurement_plan_polarizer_CN_2021_final.pdf : T. Higuchi's summary of measurement plans prepared before the beamtime
   
- notes/ : miscellaneous notes for analysis
	- notes/analysis_lod.md : the analysis log. fill this in when you make major work of analysis
	- notes/analysis_memo.md : use this to make notes of useful information for the analysis
	- notes/commands.md : list of comamnds useful for analysis 
	- notes/run_list.md : table which list information of all runs (to be filled as the analysis progresses)
	- notes/scan_list.md : table which list information of all scans (to be filled as the analysis progresses)

- refs/ : place for relevant references

- results/ : use this to place your analysis results

- tools/ : a place for codes and C++ headers needed for analysis. Edit the codes in this directory with care  
	- tools/210713_FeSi : codes which were originally placed in /home/nop/data/210713_FeSi of niki
	- tools/bin : copied from /home/nop/bin/ of niki
	- tools/ichikawa : codes provided by G. Ichikawa which can be used to split a .root file of scans using sacn logs 

## How to get data 
The contents of data/ is found in the following location:

http://www.rcnp.osaka-u.ac.jp/~thiguchi/shared/data/jparc_202107/data_210721.tar.gz

Clone this repository to your computer, move data_210721.tar.gz inside polarization-foils-bl05-202107, and then extract the tar file.



