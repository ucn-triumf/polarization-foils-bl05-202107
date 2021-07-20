#!/bin/bash
ROOTFILE_NUM="20210714000204"
SCANFILE_NAME="scan20210713_x_m2_scan_1"
GATENETLOG_NAME="../gatenetlog/gatenetlog_${SCANFILE_NAME}.txt"
echo "root MakeElist.C+(\"$ROOTFILE_NUM\",\"$GATENETLOG_NAME\")"
#root 'MakeElist.C+("'"$ROOTFILE_NUM"'","'"$GATENETLOG_NAME"'")'

I_SCAN=5
echo "root DrawOneScan.C+(\"$ROOTFILE_NUM\",$I_SCAN)"
root 'DrawOneScan.C+("'"$ROOTFILE_NUM"'",'"$I_SCAN"')'
