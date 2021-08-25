import numpy as np
import pandas as pd
import os
from glob  import glob
path = "codes/cut_scans/gatenetlog_edit/"
# fname = "gatenetlog_scan20210713_flipper_agilent_scan_1.txt"
fname = "gatenetlog_scan20210713_flipper_agilent_scan_fine_5.txt"
df_log = pd.read_csv(path + fname, names=["state", "kp"], sep="\s+")
kp_start = np.array(df_log[(df_log.index>=1)&(df_log.state=='start')].kp)
kp_end = np.array(df_log[(df_log.index>=1)&(df_log.state=='end')].kp)
print (df_log[df_log.index>0])
print (kp_start.size)
print (kp_end.size)
print (kp_end- kp_start)

