import numpy as np 
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy.polynomial.polynomial as poly
import matplotlib.ticker as plticker
mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16 

color_lst = ["lightgreen","red","lightblue","orange","orchid"]
mark = ["o","s","^","v","D"]

col_idx = 1
cv_res_file = ["HT", "MX-UT-0.1",  "MX-UT-0.5"]
for i in cv_res_file:
    file = str(i) +"-neg.ods"
    df = pd.read_excel(file, engine="odf")
    scan = df.loc[:,'scanrate'].to_numpy()
    jpa = df.loc[:,'Jpa'].to_numpy()
    jpc = df.loc[:,'Jpc'].to_numpy()
    plt.plot(scan,jpa/jpc,'o-',color=color_lst[col_idx])
    col_idx += 1