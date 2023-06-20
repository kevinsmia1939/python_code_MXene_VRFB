import numpy as np 
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from scipy.signal import savgol_filter
from scipy import linalg
from scipy.sparse.linalg import splu
from scipy.sparse import identity
import scipy.sparse as sparse
from scipy.sparse import spdiags
from impedance import preprocessing
from impedance.models.circuits import Randles, CustomCircuit
from impedance.visualization import plot_nyquist
import matplotlib as mpl
import sys
import matplotlib.ticker as plticker

mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16

def search_string_in_file(file_name, string_to_search):
    line_number = 0
    list_of_results = []
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            line_number += 1
            if string_to_search in line:
                list_of_results.append((line_number, line.rstrip()))
    return list_of_results

def CV_file2df(CV_file):
    if CV_file.lower().endswith(".csv"):
        df_CV = pd.read_csv(CV_file,usecols=[0,1])
        df_CV = np.array(df_CV)
    elif CV_file.lower().endswith(".txt"):
        df_CV = pd.read_table(CV_file, sep='\t', header=None, usecols=[0,1])
        df_CV = np.array(df_CV)
    elif CV_file.lower().endswith(".par"):
        # Search for line match beginning and end of CV data and give ln number
        start_segment = search_string_in_file(CV_file, 'Definition=Segment')[0][0]
        end_segment = search_string_in_file(CV_file, '</Segment')[0][0]
        # Count file total line number
        with open(CV_file) as f:
            ln_count = sum(1 for _ in f)
        footer = ln_count-end_segment
        df_CV = pd.read_csv(CV_file,skiprows=start_segment, skipfooter=footer,usecols=[2,3],engine='python')
        df_CV = np.array(df_CV)
    else:
        raise Exception("Unknown file type, please choose .csv, .par")
    return df_CV

def get_CV_init(df_CV):
    cv_size = df_CV.shape[0]
    volt = df_CV[:,0]
    current = df_CV[:,1]
    return cv_size, volt, current

def get_CV_peak(inv_peak_trough,cv_size, volt, current, peak_range, peak_pos, trough_pos, jpa_lns, jpa_lne, jpc_lns, jpc_lne, peak_defl_bool, trough_defl_bool):
    # If peak range is given as 0, then peak is just where peak position is
    trough_range = peak_range
    if peak_defl_bool == 1:
        peak_range = 0
        peak_curr = current[peak_pos]
        peak_volt = volt[peak_pos]   
        low_range_peak = peak_pos
        high_range_peak = peak_pos
    # Search for peak between peak_range.     
    else:
        high_range_peak = np.where((peak_pos+peak_range)>=(cv_size-1),(cv_size-1),peak_pos+peak_range)
        low_range_peak = np.where((peak_pos-peak_range)>=0,peak_pos-peak_range,0)
        peak_curr_range = current[low_range_peak:high_range_peak]
        if inv_peak_trough == False:
            peak_curr = max(peak_curr_range)
        else:
            peak_curr = min(peak_curr_range)
        peak_idx = np.argmin(np.abs(peak_curr_range-peak_curr))     
        peak_volt = volt[low_range_peak:high_range_peak][peak_idx]
        
    if trough_defl_bool == 1:
        trough_range = 0
        trough_curr = current[trough_pos]
        trough_volt = volt[trough_pos]
        high_range_trough = trough_pos
        low_range_trough = trough_pos      
    else:    
        high_range_trough = np.where((trough_pos+trough_range)>=(cv_size-1),(cv_size-1),trough_pos+trough_range)
        low_range_trough = np.where((trough_pos-trough_range)>=0,trough_pos-trough_range,0)
        trough_curr_range = current[low_range_trough:high_range_trough]
        if inv_peak_trough == False:
            trough_curr = min(trough_curr_range)
        else:
            trough_curr = max(trough_curr_range)
        trough_idx = np.argmin(np.abs(trough_curr_range-trough_curr))
        trough_volt = volt[low_range_trough:high_range_trough][trough_idx] 
    
    # If the extrapolation coordinate overlapped, just give horizontal line
    if (volt[jpa_lns:jpa_lne]).size == 0:
        volt_jpa = np.array([0, 1])
        current_jpa = np.array([0, 0])
    else:
        volt_jpa = volt[jpa_lns:jpa_lne]
        current_jpa = current[jpa_lns:jpa_lne]
        
    if (volt[jpc_lns:jpc_lne]).size == 0:
        volt_jpc = np.array([0, 1])
        current_jpc = np.array([0, 0])
    else:
        volt_jpc = volt[jpc_lns:jpc_lne]
        current_jpc = current[jpc_lns:jpc_lne]

    jpa_lnfit_coef = np.polyfit(volt_jpa,current_jpa, 1) # 1 for linear fit
    jpc_lnfit_coef = np.polyfit(volt_jpc,current_jpc, 1)
    jpa_poly1d = np.poly1d(jpa_lnfit_coef)
    jpc_poly1d = np.poly1d(jpc_lnfit_coef)
    jpa = peak_curr - jpa_poly1d(peak_volt)
    jpc = jpc_poly1d(trough_volt) - trough_curr
    return low_range_peak, high_range_peak, peak_volt, peak_curr, low_range_trough, high_range_trough, trough_volt, trough_curr, jpa, jpc, jpa_poly1d, jpc_poly1d

def cv_peak_sqrt(title,a,b, file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos):
    for i in scan_list:
        file_name = file_name1+str(i)+file_name2
        df_CV = CV_file2df(file_name)
        cv_size , volt, current = get_CV_init(df_CV)
        cy_e = volt.size
        cy_s = cy_e - 2000
        volt = volt[cy_s:cy_e]
        current = current[cy_s:cy_e]
        volt = volt[~np.isnan(volt)]
        current = current[~np.isnan(current)]        
        current = current/5*1000

        ax[a,b].plot(volt,current,'-',linewidth=0.9,label=str(i)+' mV/s')
        ax[a,b].set_xlabel('Voltage (V) vs Ag/AgCl')
        ax[a,b].set_ylabel('Current density (mA/cm$^2$)')

        
        if a == 0:
            ax[a,b].set_xlim(-0.75,0)
        if a == 1:
            ax[a,b].set_xlim(0,1.2)
        if a == 0 and b == 0:
            ax[a,b].legend(fontsize=9)
        if a == 1 and b == 0:
            ax[a,b].legend(fontsize=9)
        # elif b == 1:
        #     ax[a,b].set_xlim(0,1.2)
        # elif a == 1:
        #     ax[a,b].set_xlim(0,1.2)
        #     if b == 0:
        #         ax[a,b].legend(fontsize=9)
        # else:
        #     ax[a,b].set_xlim(0,1.2)
        # ax[a,b].set_title(title,y=1,x=0.05,pad=-14)
        
        # if b==1 and a==1:
        #     ax[a,b].set_yticks([-10,-8,-6,-4,-2,0,2,4,6,8,10,12])
        loc = plticker.MultipleLocator(base=2.0) # this locator puts ticks at regular intervals
        ax[a,b].yaxis.set_major_locator(loc)
    return

# mpl.rcParams['figure.dpi'] = 150

fig, ax = plt.subplots(2, 4, figsize=(20,10))


file_name1 = 'CV and EIS/ut/ut_0to-0.75_'
file_name2 = 'mvs.par'
scan_list = [10,9,8,7,6,5,4,3,2]
elec_area = 5
jpa_lns = 1050
jpa_lne = 1100
jpc_lns = 450
jpc_lne = 500
peak_range = 100
peak_pos = 1400
trough_pos = 900
a,b = 0,0
title = 'a)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# print(peak_sep_list1)

file_name1 = 'CV and EIS/wetted-ace-UT/UT-'
file_name2 = 'mvs-0to1.2.par'
scan_list = [50,30,20,10,5,3]
elec_area = 5
jpa_lns = 400
jpa_lne = 450
jpc_lns = 1310
jpc_lne = 1360
peak_range = 200
peak_pos = 800
trough_pos = 1600
a,b = 1,0
title = 'b)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# print(trough_volt2)
# print(peak_sep_list2)

file_name1 = 'CV and EIS/MX-HT_morescanrates/HT/HT-comp/cvmc_HT_0to-0.75_'
file_name2 = 'mvs.par'
scan_list = [10,9,8,7,6,5,4,3,2]
elec_area = 5
jpa_lns = 1680
jpa_lne = 1800
jpc_lns = 300
jpc_lne = 400
peak_range = 200
peak_pos = 1600
trough_pos = 900
a,b = 0,1
title = 'c)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)

# print(peak_sep_list3)

file_name1 = 'CV and EIS/MX-HT-1mgcm/mc3_cv_HT_1mgcm2_0to1.2_'
file_name2 = 'mvs.par'
scan_list = [50,30,20,10,5,3]
elec_area = 5
jpa_lns = 400
jpa_lne = 450
jpc_lns = 1310
jpc_lne = 1350
peak_range = 150
peak_pos = 670
trough_pos = 1600
a,b = 1,1
title = 'd)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# print(trough_volt4)
# print(peak_sep_list4)

file_name1 = 'CV and EIS/MX-UT0.1/mx-ut0.1-0to-0.75-'
file_name2 = 'mvs.par'
scan_list = [10,9,8,7,6,5,4,3,2]
elec_area = 5
jpa_lns = 1600
jpa_lne = 1700
jpc_lns = 500
jpc_lne = 600
peak_range = 200
peak_pos = 1400
trough_pos = 740
a,b = 0,2
title = 'e)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)

file_name1 = 'CV and EIS/MX-UT0.1/mx-ut0.1-0to-1.2-'
file_name2 = 'mvs.par'
scan_list = [50,30,20,10,5,3]
elec_area = 5
jpa_lns = 400
jpa_lne = 450
jpc_lns = 1310
jpc_lne = 1330
peak_range = 200
peak_pos = 900
trough_pos = 1600
a,b = 1,2
title = 'f)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)

file_name1 = 'CV and EIS/mx-ut-0.5/mx-0.5/mx-ut-0.5_0to-0.75_'
file_name2 = 'mvs.par'
scan_list = [10,9,8,7,6,5,4,3,2]
elec_area = 5
jpa_lns = 1600
jpa_lne = 1700
jpc_lns = 500
jpc_lne = 600
peak_range = 200
peak_pos = 1400
trough_pos = 740
a,b = 0,3
title = 'g)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)

file_name1 = 'CV and EIS/mx-ut-0.5/mx-0.5/mx-ut-0.5_0to1.2_'
file_name2 = 'mvs.par'
scan_list = [50,30,20,10,5,3]
elec_area = 5
jpa_lns = 400
jpa_lne = 450
jpc_lns = 1310
jpc_lne = 1330
peak_range = 200
peak_pos = 900
trough_pos = 1600
a,b = 1,3
title = 'h)'
cv_peak_sqrt(title,a,b,file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
