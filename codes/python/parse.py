import numpy as numpy
import pandas as pd
import os
from glob  import glob
os.chdir("../../data/210713_SiFe")
rootfiles = sorted(glob("20*.root"))

# print (rootfiles)
df = pd.DataFrame([])
df['fname'] = rootfiles
runnames = []
for i in range(len(df.index)):
    runnames.append(df['fname'][i][:14])
df['runnames'] = runnames
os.chdir("../../notes")
df.to_csv("run_list.csv")
