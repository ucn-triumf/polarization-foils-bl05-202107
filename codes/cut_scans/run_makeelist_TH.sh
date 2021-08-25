#!/bin/bash
ANA_HOME=$PWD
BIN_PATH=$ANA_HOME"/tools/bin/"
FESI_PATH=$ANA_HOME"/tools/210713_FeSi/"
CODE_PATH=$ANA_HOME"/codes/cut_scans/"
DATA_PATH=$ANA_HOME"/data/210713_SiFe/"

# Inputs for 20210717030252, Fe 90 nm scan
# ROOTFILE_NUM="20210717030252"
# SCANFILE_NAME="gatenetlog_scan20210713_flipper_agilent_scan_fine_3.txt"

# # Inputs for 20210717061701, Fe 50 nm scan
# ROOTFILE_NUM="20210717061701"
# SCANFILE_NAME="gatenetlog_scan20210713_flipper_agilent_scan_fine_5_mod.txt"

# Inputs for 20210716220736, Fe 30 nm scan
ROOTFILE_NUM="20210716220736"
SCANFILE_NAME="gatenetlog_scan20210713_flipper_agilent_scan_long_1.txt"

GATENETLOG_NAME=$ANA_HOME"/codes/cut_scans/gatenetlog_edit/${SCANFILE_NAME}"

cd $DATA_PATH
pwd
# echo "root \"$FESI_PATH\"MakeElist.C+(\"$ROOTFILE_NUM\",\"$GATENETLOG_NAME\")"
# root ''"$FESI_PATH"'MakeElist.C+("'"$ROOTFILE_NUM"'","'"$GATENETLOG_NAME"'")'

I_SCAN=3 # until 16
echo "root \"$CODE_PATH\"DrawOneScan.C+(\"$ROOTFILE_NUM\",\"$I_SCAN\")"
root -l ''"$CODE_PATH"'DrawOneScan.C+("'"$ROOTFILE_NUM"'",'"$I_SCAN"')'

