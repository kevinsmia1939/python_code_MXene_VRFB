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
from scipy.signal import savgol_filter
import statsmodels.api as sm

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
        jpa_arr = np.concatenate((jpa_arr, [jpa]))
        jpc_arr = np.concatenate((jpc_arr, [jpc]))
        peak_sep = peak_volt - trough_volt
        peak_sep_list.append(peak_sep)
        plt.plot(volt,current,'-',label=str(i)+' mV/s')
        
        plt.plot(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne],'r-',linewidth=4)
        plt.plot(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne],'r-',linewidth=4)        
        plt.plot((volt[jpa_lns],peak_volt),(current[jpa_lns],jpa_base),'--')
        plt.plot((volt[jpc_lns],trough_volt),(current[jpc_lns],jpc_base),'--')
        plt.plot((volt[low_range_peak],volt[high_range_peak]),(current[low_range_peak],current[high_range_peak]),"|", markersize = 10)
        plt.plot((volt[low_range_trough],volt[high_range_trough]),(current[low_range_trough],current[high_range_trough]),"|", markersize = 10)
        plt.annotate(text='', xy=(trough_volt,jpc_base), xytext=(trough_volt,trough_curr), arrowprops=dict(arrowstyle='<->'))
        plt.annotate(text='', xy=(peak_volt,jpa_base), xytext=(peak_volt,peak_curr), arrowprops=dict(arrowstyle='<-'))
        
        
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


def cv_peak_sqrt_blank(file_name1,file_name2,blank_name1,blank_name2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos):     
    volt_list=[]
    current_list=[]
    jpa_arr = []
    jpc_arr = []
    peak_sep_list=[]
    for i in scan_list:
        file_name = file_name1+str(i)+file_name2
        blank_name = blank_name1+str(i)+blank_name2
        df_CV = CV_file2df(file_name)
        _, volt, current = get_CV_init(df_CV)
        
        df_CV_blank = CV_file2df(blank_name)
        _, volt_blank, current_blank = get_CV_init(df_CV_blank)

        
        cy_e = volt.size
        cy_s = cy_e - 2000
        volt = volt[cy_s:cy_e]
        current = current[cy_s:cy_e]
        
        cy_e_blank = volt_blank.size
        cy_s_blank = cy_e_blank - 2000
        volt_blank = volt_blank[cy_s_blank:cy_e_blank]
        current_blank = current_blank[cy_s_blank:cy_e_blank]
        
        # volt = volt-volt_blank
        # current = current-current_blank
        # volt = volt - current*0
        
        volt = volt[~np.isnan(volt)]
        current = current[~np.isnan(current)]
        
        # volt_blank = volt_blank[~np.isnan(volt_blank)]
        # current_blank = current_blank[~np.isnan(current_blank)]

        current = current/elec_area*1000
        # current_blank = current_blank/elec_area*1000
        
        
        # volt_sub = np.array(volt) - np.array(volt_blank)
        # current_sub = current - current_blank
        plt.plot(volt,current,'-',label=str(i)+' mV/s')
        # plt.plot(volt_blank,current_blank,'--',label='Background ' + str(i)+' mV/s')
        # plt.plot(volt,current_sub,':',label='Subtracted CV')
        
        # current = current_sub
        low_range_peak, high_range_peak, peak_volt, peak_curr, low_range_trough, high_range_trough, trough_volt, trough_curr = get_CV_peak(volt,current,peak_range, peak_pos, trough_pos)
        
        jpa_lnfit = np.polyfit(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne], 1)
        jpa_base = jpa_lnfit[0]*peak_volt + jpa_lnfit[1]
        jpc_lnfit = np.polyfit(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne], 1)
        jpc_base = jpc_lnfit[0]*trough_volt + jpc_lnfit[1]        
        jpa = peak_curr - jpa_base
        jpc = jpc_base - trough_curr
        jpa_arr = np.concatenate((jpa_arr, [jpa]))
        jpc_arr = np.concatenate((jpc_arr, [jpc]))
        peak_sep = peak_volt - trough_volt
        peak_sep_list.append(peak_sep)


        # plt.plot(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne],'r-',linewidth=4)
        # plt.plot(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne],'r-',linewidth=4)        
        # plt.plot((volt[jpa_lns],peak_volt),(current[jpa_lns],jpa_base),'--')
        # plt.plot((volt[jpc_lns],trough_volt),(current[jpc_lns],jpc_base),'--')
        # plt.plot((volt[low_range_peak],volt[high_range_peak]),(current[low_range_peak],current[high_range_peak]),"|", markersize = 10)
        # plt.plot((volt[low_range_trough],volt[high_range_trough]),(current[low_range_trough],current[high_range_trough]),"|", markersize = 10)
        # plt.annotate(text='', xy=(trough_volt,jpc_base), xytext=(trough_volt,trough_curr), arrowprops=dict(arrowstyle='<->'))
        # plt.annotate(text='', xy=(peak_volt,jpa_base), xytext=(peak_volt,peak_curr), arrowprops=dict(arrowstyle='<-'))
        
        
        top_volt = volt[np.argmin(volt):-1]
        top_curr = current[np.argmin(volt):-1]
        
        df = pd.DataFrame(data={'volt': top_volt, 'current': top_curr})
        df = df.sort_values(['volt', 'current'], ascending=[True, False])
        der = np.gradient(df['current'],df['volt'])/30
        lowess = sm.nonparametric.lowess(der, top_volt, frac=0.1)
        smh_volt = lowess[:, 0]
        smh_curr = lowess[:, 1]
        plt.plot(top_volt,top_curr)
        plt.plot(top_volt,der,'--')
        
        plt.plot(smh_volt, smh_curr)
        
        der2 = np.gradient(smh_curr,smh_volt)/30
        # plt.plot(smh_volt,der2)
        
        lowess2 = sm.nonparametric.lowess(der2, smh_volt, frac=0.05)
        smh_volt2 = lowess2[:, 0]
        smh_curr2 = lowess2[:, 1]
        plt.plot(smh_volt2, smh_curr2) 
        plt.grid()
        
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


# fig, ax = plt.subplots(figsize=(10,10))
# file_name1_1 = 'CV and EIS/MX-HT_morescanrates/MX0.1mg-UT-0to-0.75_'
# file_name2_1 = 'mvs.par'
# blank_name1_1 = 'CV and EIS/MX-HT_morescanrates/blank/mx-ut-0.1_0to-0.75_'
# blank_name2_1 = 'mvs.par'
# scan_list1 = [10,9,8,7,6,5,4,3,2,1]
# # scan_list1 = [3]
# elec_area = 5
# jpa_lns = 1050
# jpa_lne = 1100
# jpc_lns = 480
# jpc_lne = 550
# peak_range = 100
# peak_pos = 1400
# trough_pos = 900
# volt_arr1,current_arr1,scan_rate_sqrt1,jpa_arr1,jpc_arr1,lin_fit_jpa1,lin_fit_jpc1,peak_volt1,trough_volt1,jpa_base1,jpc_base1,peak_sep_list1 = cv_peak_sqrt_blank(file_name1_1,file_name2_1, blank_name1_1, blank_name2_1, scan_list1,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.title('MX-UT-0.1mg/cm2', fontsize=16)
# plt.legend(fontsize=15)
# plt.show()
# plt.cla()
# plt.clf()

fig, ax = plt.subplots(figsize=(10,10))
file_name1_1 = 'CV and EIS/MX-HT_morescanrates/MX1mg-UT-0to-0.75_'
file_name2_1 = 'mvs.par'
scan_list1 = [10,9,8,7,6,5,4,3]
# scan_list = [10,5,3]
elec_area = 5
jpa_lns = 1050
jpa_lne = 1100
jpc_lns = 490
jpc_lne = 550
peak_range = 100
peak_pos = 1400
trough_pos = 900
volt_arr1,current_arr1,scan_rate_sqrt1,jpa_arr1,jpc_arr1,lin_fit_jpa1,lin_fit_jpc1,peak_volt1,trough_volt1,jpa_base1,jpc_base1,peak_sep_list1 = cv_peak_sqrt(file_name1_1,file_name2_1,scan_list1,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
plt.title('MX-UT-0.1mg/cm2', fontsize=16)
plt.legend(fontsize=15)
plt.show()
plt.cla()
plt.clf()

# fig, ax = plt.subplots(figsize=(10,10))
# file_name1_2 = 'CV and EIS/MX-HT_morescanrates/MX0.05mg-UT-0to-0.75_'
# file_name2_2 = 'mvs.par'
# scan_list2 = [10,8,7,6,5,4,3]
# # scan_list2 = [10,5,3]
# elec_area = 5
# jpa_lns = 1155
# jpa_lne = 1160
# jpc_lns = 500
# jpc_lne = 560
# peak_range = 100
# peak_pos = 1400
# trough_pos = 850
# volt_arr2,current_arr2,scan_rate_sqrt2,jpa_arr2,jpc_arr2,lin_fit_jpa2,lin_fit_jpc2,peak_volt2,trough_volt2,jpa_base2,jpc_base2,peak_sep_list2 = cv_peak_sqrt(file_name1_2,file_name2_2,scan_list2,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.title('MX-HT-1mg/cm2', fontsize=16)
# plt.legend(fontsize=15)
# plt.show()
# plt.cla()
# plt.clf()

fig, ax = plt.subplots(figsize=(25,25))
file_name1_3 = 'CV and EIS/MX-HT_morescanrates/HT/HT-comp/cvmc_HT_0to-0.75_'
file_name2_3 = 'mvs.par'
blank_name1_3 = 'CV and EIS/MX-HT_morescanrates/blank/ht_0to-0.75_'
blank_name2_3 = 'mvs.par'
# scan_list = [10,9,8,7,6,5,4,3,2]
# scan_list = [4,3,2]
scan_list = [4]
elec_area = 5
jpa_lns = 1050
jpa_lne = 1100
jpc_lns = 450
jpc_lne = 500
peak_range = 100
peak_pos = 1400
trough_pos = 900
volt_arr3,current_arr3,scan_rate_sqrt3,jpa_arr3,jpc_arr3,lin_fit_jpa3,lin_fit_jpc3,peak_volt3,trough_volt3,jpa_base3,jpc_base3,peak_sep_list3 = cv_peak_sqrt_blank(file_name1_3,file_name2_3,blank_name1_3,blank_name2_3,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.title('HT', fontsize=16)
plt.legend(fontsize=15)
plt.show()
plt.cla()
plt.clf()

# fig, ax = plt.subplots(figsize=(10,10))
# file_name1_5 = 'CV and EIS/MX-HT_morescanrates/MX0.0mg-UT-0to-0.75_'
# file_name2_5 = 'mvs.par'
# blank_name1_5 = 'CV and EIS/MX-HT_morescanrates/blank/ht_0to-0.75_'
# blank_name2_5 = 'mvs.par'
# scan_list = [10,9,8,7,6,5,4,3,2,1]
# # scan_list = [4,3,2,1]
# # scan_list = [4]
# elec_area = 5
# jpa_lns = 1050
# jpa_lne = 1100
# jpc_lns = 450
# jpc_lne = 500
# peak_range = 100
# peak_pos = 1400
# trough_pos = 900
# volt_arr5,current_arr5,scan_rate_sqrt5,jpa_arr5,jpc_arr5,lin_fit_jpa5,lin_fit_jpc5,peak_volt5,trough_volt5,jpa_base5,jpc_base5,peak_sep_list5 = cv_peak_sqrt_blank(file_name1_5,file_name2_5,blank_name1_5,blank_name2_5,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# # plt.title('HT', fontsize=16)
# plt.legend(fontsize=15)
# plt.show()
# plt.cla()
# plt.clf()


# fig, ax = plt.subplots(figsize=(10,10))



# plt.plot(np.sqrt(np.flip(scan_list1)),jpa_arr1,'bo-',label='MX-UT')
# plt.plot(np.sqrt(np.flip(scan_list1)),-jpc_arr1,'bo-')


# jpa_arr2 = np.array([0.2,0.6])
# # print(np.sqrt(np.flip(scan_list))[0:2])

# plt.plot(np.sqrt(np.flip(scan_list))[0:2],jpa_arr2,'rs-',label='MX-HT')
# # plt.plot(np.sqrt(np.flip(scan_list))[0:2],-jpc_arr2[0:2],'rs-')

# jpa_arr3 = [0.2,0.4]
# plt.plot(np.sqrt(np.flip(scan_list))[0:2],jpa_arr3,'g^-',label='HT')
# plt.plot(np.sqrt(np.flip(scan_list))[0:2],-jpc_arr3[0:2],'g^-')

# plt.legend(fontsize=15)
# plt.xlabel('Scan rate$^{1/2}$ (mV/s)$^{1/2}$')
# plt.ylabel('Peak current (mA)')
# plt.show()
# plt.cla()
# plt.clf()

# # Reversibility
# fig, ax = plt.subplots(figsize=(10,10))
# plt.plot(np.flip(scan_list1),jpa_arr1/jpc_arr1,'bo-',label='MX-UT')
# # plt.plot(np.flip(scan_list)[0:2],(jpa_arr2[0:2]/jpc_arr2[0:2].T).T,'rs-',label='MX-HT')
# plt.plot(np.flip(scan_list)[0:2],(jpa_arr3[0:2]/jpc_arr3[0:2].T).T,'g^-',label='HT')
# plt.legend(fontsize=15)
# plt.xlabel('Scan rate (mV/s)')
# plt.ylabel('jpa/jpc')
# plt.show()
# plt.cla()
# plt.clf()

# fig, ax = plt.subplots(figsize=(10,10))
# plt.plot(scan_list1,peak_sep_list1,'o-',label='MX-UT')
# # plt.plot(scan_list[1:3],peak_sep_list2[1:3],'s-',label='MX-HT')
# # plt.plot(peak_sep_list3[1:3],scan_list[1:3],'^-',label='HT')
# ht_x_3 = -0.4247881 - -0.6707667
# ht_x_5 = -0.4165071 - -0.6913159
# plt.plot(scan_list[1:3],[ht_x_5,ht_x_3],'^-',label='HT')
# plt.legend(fontsize=15)
# plt.ylabel('Peak separation (V)')
# plt.xlabel('Scan rate (mV/s)')
# plt.show()
