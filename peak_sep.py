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

fig, ax = plt.subplots(figsize=(15,15))


volt_list=[]
current_list=[]
for i in [20,10,5,3]:
    file_name = 'CV and EIS/metal_Ti/cv_Ti_0to1.2_'+str(i)+'mvs.par'
    df_CV = CV_file2df(file_name)
    _, volt, current = get_CV_init(df_CV)
    volt = volt - current*1.45
    current = current/1.2*1000
    volt_list.append(volt)
    current_list.append(current)
    
current_arr = np.array(current_list)
volt_arr = np.array(volt_list)

# plt.plot(volt_arr,current_arr)

jpa_lns = 100
jpa_lne = 200
jpc_lns = 1210
jpc_lne = 1280
jpa_list = []
jpc_list = []
# trough_curr_list=[]
# peak_curr_list=[]
peak_sep_list=[]
for cycle in np.arange(0,4,1):
    low_range_peak, high_range_peak, peak_volt, peak_curr, low_range_trough, high_range_trough, trough_volt, trough_curr = get_CV_peak(volt_arr[cycle],current_arr[cycle],190,820,1600)
    
    jpa_lnfit = np.polyfit(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne], 1)
    jpa_base = jpa_lnfit[0]*peak_volt + jpa_lnfit[1]

    jpc_lnfit = np.polyfit(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne], 1)
    jpc_base = jpc_lnfit[0]*trough_volt + jpc_lnfit[1]
    
    jpa = peak_curr - jpa_base
    jpc = jpc_base - trough_curr
    
    jpa_list.append(jpa)
    jpc_list.append(jpc)
    # trough_curr_list.append(trough_volt)
    # peak_curr_list.append(peak_volt)
    peak_sep = peak_volt - trough_volt
    peak_sep_list.append(peak_sep)
    # plt.plot((volt_single[cycle][jpa_lns],peak_volt),(current_single[cycle][jpa_lns],jpa_base),'--')
    # plt.plot((volt_single[cycle][jpc_lns],trough_volt),(current_single[cycle][jpc_lns],jpc_base),'--')
    
    # plt.plot(peak_volt,peak_curr,'o')
    # plt.plot(trough_volt,trough_curr,'s')
    
    plt.annotate(text='', xy=(trough_volt,jpc_base), xytext=(trough_volt,trough_curr), arrowprops=dict(arrowstyle='<->'))
    plt.annotate(text='', xy=(peak_volt,jpa_base), xytext=(peak_volt,peak_curr), arrowprops=dict(arrowstyle='<-'))
    
    # plt.plot((volt_arr[cycle][low_range_peak],volt_arr[cycle][high_range_peak]),(current_arr[cycle][low_range_peak],current_arr[cycle][high_range_peak]),"|", markersize = 10)
    # plt.plot((volt_arr[cycle][low_range_trough],volt_arr[cycle][high_range_trough]),(current_arr[cycle][low_range_trough],current_arr[cycle][high_range_trough]),"|", markersize = 10)
    plt.plot(volt_arr[cycle],current_arr[cycle])
    
plt.xlim(0,1.2)
plt.ylim(-10,10)
plt.show()

jpa_list = np.flip(np.vstack(jpa_list))
jpc_list = np.flip(np.vstack(jpc_list))
peak_sep_list = np.flip(np.vstack(peak_sep_list))
# jpa_list.reverse()
# jpc_list.reverse()
# trough_curr_list.reverse()
# peak_curr_list.reverse()
# peak_sep_list.reverse()

fig, ax = plt.subplots(figsize=(8,8))
number_cycle = [3,5,10,20]
plt.plot(number_cycle,np.array(jpa_list),'^',label='$V^{4+}$ → $V^{5+}$')
plt.plot(number_cycle,-np.array(jpc_list),'o',label='$V^{5+}$ → $V^{4+}$')
plt.xlabel('Scan rate', fontsize=16)
plt.ylabel('Peak current density $(mA/cm^2)$', fontsize=16)
plt.legend(fontsize=15)
plt.title('UT', fontsize=16)
plt.show()

fig, ax = plt.subplots(figsize=(8,8))
plt.plot(number_cycle,np.array(peak_sep_list)*1000,'o')
# plt.plot(number_cycle,-np.array(trough_curr_list),'o',label='$V^{5+}$ → $V^{4+}$')
plt.xlabel('Scan rate', fontsize=16)
plt.ylabel('Peak separation $(mV$)', fontsize=16)
# plt.legend(fontsize=15)
plt.title('UT', fontsize=16)
plt.show()

# df_CV1 = CV_file2df('CV and EIS/metal_Ti/cv_Ti_0to1.2_50mvs.par')
# _, volt1, current1 = get_CV_init(df_CV1)
# plt.plot(volt1, current1/1.2*1000,label='50 mV/s')
# low_range_peak, high_range_peak, peak_volt, peak_curr, low_range_trough, high_range_trough, trough_volt, trough_curr = get_CV_peak(volt1,current1,100,650,1500)


# df_CV2 = CV_file2df('CV and EIS/metal_Ti/cv_Ti_0to1.2_30mvs.par')
# _, volt2, current2 = get_CV_init(df_CV2)
# plt.plot(volt2, current2/1.2*1000,label='30 mV/s')

# df_CV3 = CV_file2df('CV and EIS/metal_Ti/cv_Ti_0to1.2_20mvs.par')
# _, volt3, current3 = get_CV_init(df_CV3)
# plt.plot(volt3, current3/1.2*1000,label='20 mV/s')

# df_CV4 = CV_file2df('CV and EIS/metal_Ti/cv_Ti_0to1.2_10mvs.par')
# _, volt4, current4 = get_CV_init(df_CV4)
# plt.plot(volt4, current4/1.2*1000,label='10 mV/s')

# df_CV5 = CV_file2df('CV and EIS/metal_Ti/cv_Ti_0to1.2_5mvs.par')
# _, volt5, current5 = get_CV_init(df_CV5)
# plt.plot(volt5, current5/1.2*1000,label='5 mV/s')

# df_CV6 = CV_file2df('CV and EIS/metal_Ti/cv_Ti_0to1.2_3mvs.par')
# _, volt6, current6 = get_CV_init(df_CV6)
# plt.plot(volt6, current6/1.2*1000,label='3 mV/s')

# plt.xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
# plt.ylabel('Current density (mA/cm$^2$)', fontsize=16)
# # plt.legend(fontsize=15)
# plt.xlim(0,1.2)
# # plt.ylim(-10,10)
# plt.show()
# plt.cla()
# plt.clf()

