import numpy as numpy
import pandas as pd
import os
from glob  import glob
from pathlib import Path

#os.chdir("../../data/202107013_Polarizer")
#scanfiles = glob("gatenetlog_scan20210713_*")
os.chdir("../../data/202107013_Polarizer")
scanfiles = glob("gatenetlog_scan20210713_*")

scanfiles_sorted = sorted(scanfiles, key=os.path.getmtime)
# print (rootfiles)
df = pd.DataFrame([])
df['fname'] = scanfiles_sorted 
scannames = []
for i in range(len(df.index)):
    scannames.append(df['fname'][i][11:])
df['scanname'] = scannames
os.chdir("../../notes")
df.to_csv("scan_list.csv")
