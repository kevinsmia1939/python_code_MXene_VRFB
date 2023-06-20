import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16 
fig, ax = plt.subplots(figsize=(6,12))

xps_lst = ['UT','HT','MX-UT-0.1','MX-UT-0.5']

shift_idx = 0
shift=[0,-20000,70000,180000]
shift_l = 0
for xps in xps_lst:
    print(shift[shift_idx])
    excel = '../Excel sheets/' + xps + '.xlsx'
    df = pd.read_excel(excel,sheet_name='0-wide')  
    bind = df.iloc[1:,0].to_numpy() #Put C1s to 284.5 eV
    count = df.iloc[1:,1].to_numpy()-shift_l+shift[shift_idx]
    plt.plot(bind,count,linewidth = 1,color='black')
    plt.text(1000, -shift_l+shift[shift_idx]+300000, xps, rotation=0, fontsize=18,horizontalalignment='center',verticalalignment='center')
    shift_l += 460000
    shift_idx += 1

plt.xlabel('Binding energy (eV)')
plt.ylabel('Intensity (arb. unit)')
plt.xlim(bind[0],bind[-1])
plt.text(284.5, 460000, 'C1s', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(450, -680000, 'Ti2p', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(450, -1000000, 'Ti2p', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(531, 140000, 'O1s', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(531, -340000, 'O1s', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(531, -700000, 'O1s', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(531, -1040000, 'O1s', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(200, -780000, 'Cl2p', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(200, -1130000, 'Cl2p', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(685, -710000, 'F1s', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.text(685, -1040000, 'F1s', rotation=90, fontsize=12,horizontalalignment='center',verticalalignment='center')
plt.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)

