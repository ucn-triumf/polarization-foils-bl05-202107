#!/bin/bash
ANA_HOME=$PWD
BIN_PATH=$ANA_HOME"/tools/bin/"
FESI_PATH=$ANA_HOME"/tools/210713_FeSi/"
CODE_PATH=$ANA_HOME"/codes/cut_scans/"
DATA_PATH=$ANA_HOME"/data/210713_SiFe/"

# Inputs for 20210716210153, Fe 30 nm scan (1)
# ROOTFILE_NUM="20210716210153"
# SCANFILE_NAME="gatenetlog_scan20210713_flipper_agilent_scan_1.txt"

# # Inputs for 20210716220736, Fe 30 nm scan (2)
# ROOTFILE_NUM="20210716220736"
# SCANFILE_NAME="gatenetlog_scan20210713_flipper_agilent_scan_long_1.txt"


# # Inputs for 20210717030252, Fe 90 nm scan
# ROOTFILE_NUM="20210717030252"
# SCANFILE_NAME="gatenetlog_scan20210713_flipper_agilent_scan_fine_3.txt"

# # Inputs for 20210717061701, Fe 50 nm scan
ROOTFILE_NUM="20210717061701"
SCANFILE_NAME="gatenetlog_scan20210713_flipper_agilent_scan_fine_5_mod.txt"


GATENETLOG_NAME=$ANA_HOME"/codes/cut_scans/gatenetlog_edit/${SCANFILE_NAME}"
cd $DATA_PATH
pwd
# echo "root \"$FESI_PATH\"MakeElist.C+(\"$ROOTFILE_NUM\",\"$GATENETLOG_NAME\")"
# root -l ''"$FESI_PATH"'MakeElist.C+("'"$ROOTFILE_NUM"'","'"$GATENETLOG_NAME"'")'

I_SCAN=6
 ## until 9 for run 20210716210153, 
 ## until 3 for run 20210716220736
 ## until 13 for run 20210717030252
 ## until 16 for run 20210717061701. index=16 is a incomplete data

echo "root \"$CODE_PATH\"DrawOneScan.C+(\"$ROOTFILE_NUM\",\"$I_SCAN\")"
root -l ''"$CODE_PATH"'DrawOneScan.C+("'"$ROOTFILE_NUM"'",'"$I_SCAN"')'

cd $ANA_HOME

