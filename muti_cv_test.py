import numpy as np 
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

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

def get_CV_peak(volt,current,peak_range, peak_pos, trough_pos):
    # Search for peak between peak_range.
    cv_size = volt.shape[0]
    # volt, current = get_CV_init(df_CV)  

    high_range_peak = np.where((peak_pos+peak_range)>=(cv_size-1),(cv_size-1),peak_pos+peak_range)
    low_range_peak = np.where((peak_pos-peak_range)>=0,peak_pos-peak_range,0)
    peak_curr_range = current[low_range_peak:high_range_peak]
    peak_curr = max(peak_curr_range)
    peak_idx = np.argmin(np.abs(peak_curr_range-peak_curr))
    peak_volt = volt[low_range_peak:high_range_peak][peak_idx]

    high_range_trough = np.where((trough_pos+peak_range)>=(cv_size-1),(cv_size-1),trough_pos+peak_range)
    low_range_trough = np.where((trough_pos-peak_range)>=0,trough_pos-peak_range,0)
    trough_curr_range = current[low_range_trough:high_range_trough]
    trough_curr = min(trough_curr_range)
    trough_idx = np.argmin(np.abs(trough_curr_range-trough_curr))
    trough_volt = volt[low_range_trough:high_range_trough][trough_idx]  
    return low_range_peak, high_range_peak, peak_volt, peak_curr, low_range_trough, high_range_trough, trough_volt, trough_curr

def get_CV(df_CV,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_volt,trough_volt):
    # Select the points to extrapolate.
    if jpa_lns == jpa_lne:
        jpa_lne = jpa_lns+1
    if jpa_lns > jpa_lne:
        save_val_jpa = jpa_lns
        jpa_lns = jpa_lne
        jpa_lne = save_val_jpa
    if jpc_lns == jpc_lne:
        jpc_lne = jpc_lns+1
    if jpc_lns > jpc_lne:
        save_val_jpc = jpc_lns
        jpc_lns = jpc_lne
        jpc_lne = save_val_jpc
        
    cv_size, volt, current = get_CV_init(df_CV)    

    jpa_lnfit = np.polyfit(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne], 1)
    jpa_base = jpa_lnfit[0]*peak_volt + jpa_lnfit[1]

    jpc_lnfit = np.polyfit(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne], 1)
    jpc_base = jpc_lnfit[0]*trough_volt + jpc_lnfit[1]
    return jpa_lns,jpa_lne,jpc_lns,jpc_lne, volt, current, jpa_base, jpc_base

def cv_peak_sqrt(file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne):     
    volt_list=[]
    current_list=[]
    jpa_arr = []
    jpc_arr = []
    peak_sep_list=[]
    for i in scan_list:
        file_name = file_name1+str(i)+file_name2
        df_CV = CV_file2df(file_name)
        _, volt, current = get_CV_init(df_CV)
        volt = volt - current*0
        current = current/elec_area*1000
        
        low_range_peak, high_range_peak, peak_volt, peak_curr, low_range_trough, high_range_trough, trough_volt, trough_curr = get_CV_peak(volt,current,200,700,1600)
        
        # _,_,_,_,_,_, jpa_base, jpc_base = get_CV(df_CV,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_volt,trough_volt)
        
        jpa_lnfit = np.polyfit(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne], 1)
        jpa_base = jpa_lnfit[0]*peak_volt + jpa_lnfit[1]
        jpc_lnfit = np.polyfit(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne], 1)
        jpc_base = jpc_lnfit[0]*trough_volt + jpc_lnfit[1]
        
        jpa = peak_curr - jpa_base
        jpc = jpc_base - trough_curr
        
        jpa_arr = np.concatenate((jpa_arr, [jpa]))
        jpc_arr = np.concatenate((jpc_arr, [jpc]))
        # jpa_arr.append(jpa)
        # jpc_arr.append(jpc)
        peak_sep = peak_volt - trough_volt
        peak_sep_list.append(peak_sep)
        
        volt_list.append(volt)
        current_list.append(current)
  
    current_arr = np.array(current_list)
    volt_arr = np.array(volt_list)
    jpa_arr = np.flip(np.vstack(jpa_arr))
    jpc_arr = np.flip(np.vstack(jpc_arr))
    peak_sep_list = np.vstack(peak_sep_list)
    
    scan_rate = np.flip(np.array(scan_list))
    scan_rate_sqrt = scan_rate**0.5
    lin_fit_jpa = np.polyfit(scan_rate_sqrt,jpa_arr, 1)
    lin_fit_jpc = np.polyfit(scan_rate_sqrt,-jpc_arr, 1)
    return volt_arr,current_arr,scan_rate_sqrt,jpa_arr,jpc_arr,lin_fit_jpa,lin_fit_jpc


file_name1 = 'CV and EIS/metal_Ti/cv_Ti_0to1.2_'
file_name2 = 'mvs.par'
scan_list = [20,10,5,3]
elec_area = 1.2
jpa_lns = 300
jpa_lne = 400
jpc_lns = 1210
jpc_lne = 1280
volt_arr,current_arr,scan_rate_sqrt1,jpa_arr1,jpc_arr1,lin_fit_jpa1,lin_fit_jpc1 = cv_peak_sqrt(file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne)

file_name3 = 'CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_'
file_name4 = 'mVs_beaker_2.par'
scan_list = [20,10,5,3]
elec_area = 5
jpa_lns = 270
jpa_lne = 460
jpc_lns = 1280
jpc_lne = 1300
volt_arr2,current_arr2,scan_rate_sqrt2,jpa_arr2,jpc_arr2,lin_fit_jpa2,lin_fit_jpc2 = cv_peak_sqrt(file_name3,file_name4,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne)

fig, ax = plt.subplots(figsize=(10,10))
plt.plot(volt_arr2.T,current_arr2.T)
plt.show()

plt.plot(scan_rate_sqrt1,jpa_arr1,'ro',label='Ti')
plt.plot(scan_rate_sqrt1,-jpc_arr1,'ro')
plt.plot(scan_rate_sqrt1,lin_fit_jpa1[0]*scan_rate_sqrt1+lin_fit_jpa1[1],'r-')
plt.plot(scan_rate_sqrt1,lin_fit_jpc1[0]*scan_rate_sqrt1+lin_fit_jpc1[1],'r-')

plt.plot(scan_rate_sqrt2,jpa_arr2,'b^',label='HT-CP')
plt.plot(scan_rate_sqrt2,-jpc_arr2,'b^')
plt.plot(scan_rate_sqrt2,lin_fit_jpa2[0]*scan_rate_sqrt2+lin_fit_jpa2[1],'b-')
plt.plot(scan_rate_sqrt2,lin_fit_jpc2[0]*scan_rate_sqrt2+lin_fit_jpc2[1],'b-')




plt.xlabel('Scan rate$^{1/2}$ (mV/s)$^{1/2}$', fontsize=16)
plt.ylabel('Peak current density $(mA/cm^2)$', fontsize=16)
plt.legend(fontsize=15)
plt.title('UT', fontsize=16)
plt.show()
