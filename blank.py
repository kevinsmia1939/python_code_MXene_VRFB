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

def search_string_in_file(file_name, string_to_search):
    line_number = 0
    list_of_results = []
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            line_number += 1
            if string_to_search in line:
                list_of_results.append((line_number, line.rstrip()))
    return list_of_results

def whihen(y, d, lmbd):
    m = len(y)
    
    shape = (m-d, m)
    diagonals = np.zeros(2*d + 1)
    diagonals[d] = 1.
    for i in range(d):
        diff = diagonals[:-1] - diagonals[1:]
        diagonals = diff
    offsets = np.arange(d+1)    
    D = sparse.diags(diagonals, offsets, shape, format="csc")

    """
    Handle missing data points.
    Generate weight array from input data, element without value is convert to 
    0 and with value is 1.
    """
    S = np.where(np.isnan(y), y, 1)
    S = np.where(~np.isnan(S), S, 0)
    # print(S)
    W = spdiags(S, 0, m, m)

    C = W + lmbd * D.T @ D
    
    """
    scipy.splu does not able to handle NaN, so convert to 0.
    """
    y = np.where(~np.isnan(y), y, 0)
    z = splu(C).solve(y)

    return z

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
    cv_size = volt.shape[0]
    high_range_peak = np.where((peak_pos+peak_range)>=(cv_size-1),(cv_size-1),peak_pos+peak_range)
    low_range_peak = np.where((peak_pos-peak_range)>=0,peak_pos-peak_range,0)
    
    # print(low_range_peak)
    # print(high_range_peak)
    peak_curr_range = current[low_range_peak:high_range_peak]
    # print(current[low_range_peak:high_range_peak])
    # print(peak_curr_range)
    
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

def cv_peak_sqrt(file_name1,file_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos):     
    volt_list=[]
    current_list=[]
    jpa_arr = []
    jpc_arr = []
    peak_sep_list=[]
    for i in scan_list:
        file_name = file_name1+str(i)+file_name2
        df_CV = CV_file2df(file_name)
        _, volt, current = get_CV_init(df_CV)
        cy_e = volt.size
        cy_s = cy_e - 2000
        volt = volt[cy_s:cy_e]
        current = current[cy_s:cy_e]

        volt = volt[~np.isnan(volt)]
        current = current[~np.isnan(current)]
        current = current/elec_area*1000
        
        low_range_peak, high_range_peak, peak_volt, peak_curr, low_range_trough, high_range_trough, trough_volt, trough_curr = get_CV_peak(volt,current,peak_range, peak_pos, trough_pos)
        
        jpa_lnfit = np.polyfit(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne], 1)
        jpa_base = jpa_lnfit[0]*peak_volt + jpa_lnfit[1]
        jpc_lnfit = np.polyfit(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne], 1)
        jpc_base = jpc_lnfit[0]*trough_volt + jpc_lnfit[1]
        
        jpa = peak_curr - jpa_base
        jpc = jpc_base - trough_curr
        rev = jpa/jpc
        jpa_arr = np.concatenate((jpa_arr, [jpa]))
        jpc_arr = np.concatenate((jpc_arr, [jpc]))
        peak_sep = peak_volt - trough_volt
        peak_sep_list.append(peak_sep)
        
        plt.plot(volt,current,'-',label=str(i)+' mV/s')

    current_arr = 1
    volt_arr = 1
    jpa_arr = np.flip(np.vstack(jpa_arr))
    jpc_arr = np.flip(np.vstack(jpc_arr))
    peak_sep_list = np.vstack(peak_sep_list)
    
    scan_rate = np.flip(np.array(scan_list))
    scan_rate_sqrt = scan_rate**0.5
    lin_fit_jpa = np.polyfit(scan_rate_sqrt,jpa_arr, 1)
    lin_fit_jpc = np.polyfit(scan_rate_sqrt,-jpc_arr, 1)
    return volt_arr,current_arr,scan_rate_sqrt,jpa_arr,jpc_arr,lin_fit_jpa,lin_fit_jpc,peak_volt,trough_volt,jpa_base,jpc_base,peak_sep_list

def cv_diff(file_name1,file_name2,scan_list,elec_area):     
    volt_list=[]
    current_list=[]
    for i in scan_list:
        file_name = file_name1+str(i)+file_name2
        df_CV = CV_file2df(file_name)
        _, volt, current = get_CV_init(df_CV)
        cy_e = volt.size
        cy_s = cy_e - 2000
        volt = volt[cy_s:cy_e]
        current = current[cy_s:cy_e]

        volt = volt[~np.isnan(volt)]
        current = current[~np.isnan(current)]
        current = current/elec_area*1000
        

        
        vertex_min = np.argmin(volt)
        vertex_max = np.argmax(volt)
        print(vertex_min)
        print(vertex_max)
        # bot_volt = volt[0:1000]
        # bot_curr = current[0:1000]
        plt.plot(volt,current,'-',label=str(i)+' mV/s')
        
        top_volt = volt[vertex_min:1999]
        top_curr = current[vertex_min:1999]
        
        top_dydx = np.gradient(top_curr,top_volt)
        # bot_dydx = np.gradient(bot_curr,bot_volt)
        # plt.plot(x,y)


        plt.plot(top_volt,top_dydx/30,'.')
        
        smh_top_dydx = savgol_filter(top_dydx/30, 50, 3) # window size 51, polynomial order 3
        plt.plot(top_volt,smh_top_dydx, color='red')
        # plt.plot(bot_volt,bot_dydx/10,'-')
    return volt,current


# fig, ax = plt.subplots(figsize=(10,10))
# file_name1_1 = 'CV and EIS/MX-HT_morescanrates/MX0.1mg-UT-0to-0.75_'
# file_name2_1 = 'mvs.par'
# scan_list1 = [10,9,8,7,6,5,4,3]
# elec_area = 5
# jpa_lns = 1050
# jpa_lne = 1100
# jpc_lns = 490
# jpc_lne = 550
# peak_range = 100
# peak_pos = 1400
# trough_pos = 900
# volt_arr1,current_arr1,scan_rate_sqrt1,jpa_arr1,jpc_arr1,lin_fit_jpa1,lin_fit_jpc1,peak_volt1,trough_volt1,jpa_base1,jpc_base1,peak_sep_list1 = cv_peak_sqrt(file_name1_1,file_name2_1,scan_list1,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.title('MX-UT-0.1mg/cm2', fontsize=16)
# plt.legend(fontsize=15)
# plt.show()
# plt.cla()
# plt.clf()


fig, ax = plt.subplots(figsize=(15,15))
file_name1_3 = 'CV and EIS/MX-HT_morescanrates/HT/HT-comp/cvmc_HT_0to-0.75_'
file_name2_3 = 'mvs.par'
# scan_list = [10,9,8,7,6,5,4,3,2,1]
# scan_list = [4,3,2,1]
scan_list = [2]
elec_area = 5
jpa_lns = 1050
jpa_lne = 1100
jpc_lns = 450
jpc_lne = 500
peak_range = 100
peak_pos = 1400
trough_pos = 900
volt_arr3,current_arr3,scan_rate_sqrt3,jpa_arr3,jpc_arr3,lin_fit_jpa3,lin_fit_jpc3,peak_volt3,trough_volt3,jpa_base3,jpc_base3,peak_sep_list3 = cv_peak_sqrt(file_name1_3,file_name2_3,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)

volt,current = cv_diff(file_name1_3,file_name2_3,scan_list,elec_area)

plt.title('HT', fontsize=16)
plt.legend(fontsize=15)
plt.show()
plt.cla()
plt.clf()


# fig, ax = plt.subplots(figsize=(10,10))



# plt.plot(np.sqrt(np.flip(scan_list1)),jpa_arr1,'bo-',label='MX-UT')
# plt.plot(np.sqrt(np.flip(scan_list1)),-jpc_arr1,'bo-')


# jpa_arr2 = np.array([0.2,0.6])
# print(np.sqrt(np.flip(scan_list))[0:2])

# plt.plot(np.sqrt(np.flip(scan_list))[0:2],jpa_arr2,'rs-',label='MX-HT')
