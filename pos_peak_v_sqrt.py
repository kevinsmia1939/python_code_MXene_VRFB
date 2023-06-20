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
        # print(volt.shape)
        cy_e = volt.size
        cy_s = cy_e - 2000
        volt = volt[cy_s:cy_e]
        current = current[cy_s:cy_e]

        
        # volt = volt - current*0
        
        volt = volt[~np.isnan(volt)]
        current = current[~np.isnan(current)]
        
        # d=2
        # lmbd = 100000
        # volt = whihen(volt, d, lmbd)
        
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
        
        # volt_list.append(volt)
        # current_list.append(current)
        
        # plt.plot(volt[jpa_lns:jpa_lne],current[jpa_lns:jpa_lne],'r-',linewidth=4)
        # plt.plot(volt[jpc_lns:jpc_lne],current[jpc_lns:jpc_lne],'r-',linewidth=4)
        
        # plt.plot((volt[jpa_lns],peak_volt),(current[jpa_lns],jpa_base),'--')
        # plt.plot((volt[jpc_lns],trough_volt),(current[jpc_lns],jpc_base),'--')

        # plt.plot((volt[low_range_peak],volt[high_range_peak]),(current[low_range_peak],current[high_range_peak]),"|", markersize = 10)
        # plt.plot((volt[low_range_trough],volt[high_range_trough]),(current[low_range_trough],current[high_range_trough]),"|", markersize = 10)
        # plt.annotate(text='', xy=(trough_volt,jpc_base), xytext=(trough_volt,trough_curr), arrowprops=dict(arrowstyle='<->'))
        # plt.annotate(text='', xy=(peak_volt,jpa_base), xytext=(peak_volt,peak_curr), arrowprops=dict(arrowstyle='<-'))

        # vertex_min = np.argmin(volt)
        # vertex_max = np.argmax(volt)
        # print(vertex_min)
        # print(vertex_max)
        # bot_volt = volt[0:1000]
        # bot_curr = current[0:1000]
        
        # top_volt = volt[vertex_min:1999]
        # top_curr = current[vertex_min:1999]
        
        # plt.plot(top_volt,top_curr,'o')
        # print(volt.shape)
        plt.plot(volt,current,'-',label=str(i)+' mV/s')
        # yhat = savgol_filter(current, 100, 4) # window size 51, polynomial order 3
        # plt.plot(top_volt,yhat,'-')
        # plt.plot(volt,yhat,'.')
        

        # smh_top_volt = whihen(top_volt, d, lmbd)
        # smh_bot_volt = whihen(bot_volt, d, lmbd)
        # plt.plot(smh_top_volt,top_curr,'.-')

        # top_dydx = np.gradient(top_curr,smh_top_volt)
        # bot_dydx = np.gradient(bot_curr,smh_bot_volt)
        # # plt.plot(x,y)
        # plt.plot(top_volt,top_dydx/10,'-')
        # plt.plot(bot_volt,bot_dydx/10,'-')
        
        # smh_bot_dydx = whihen(bot_dydx, d, lmbd)
        
        # plt.plot(bot_volt,smh_bot_dydx/10,'-')

        # plt.ylim(-4,2)
        # plt.xlim(-0.75,-0.1)
        
        top_volt = volt[0:np.argmax(volt)]
        top_curr = current[0:np.argmax(volt)]
        df = pd.DataFrame(data={'volt': top_volt, 'current': top_curr})
        df = df.sort_values(['volt', 'current'], ascending=[True, False])
        der = np.gradient(df['current'],df['volt'])/10
        lowess = sm.nonparametric.lowess(der, top_volt, frac=0.1)
        smh_volt = lowess[:, 0]
        smh_curr = lowess[:, 1]
        plt.plot(top_volt,top_curr)
        
        # plt.plot(top_volt,der,'--')
        
        plt.plot(smh_volt, smh_curr)
        
        der2 = np.gradient(smh_curr,smh_volt)/10
        # plt.plot(smh_volt,der2)
        
        lowess2 = sm.nonparametric.lowess(der2, smh_volt, frac=0.05)
        smh_volt2 = lowess2[:, 0]
        smh_curr2 = lowess2[:, 1]
        plt.plot(smh_volt2, smh_curr2) 
        plt.grid()
        
    # print(current_list)
    # current_arr = np.array(current_list)
    # volt_arr = np.array(volt_list)
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

mpl.rcParams['figure.dpi'] = 150
# fig, ax = plt.subplots(figsize=(15,15))

# fig, ax = plt.subplots(3, 2, figsize=(20,20))


# file_name1_1 = 'CV and EIS/MX-HT-1mgcm/mc3_cv_Ti_0to1.2_'
# file_name2_1 = 'mvs.par'
# scan_list = [50,40,30,20,10,5,3]
# # scan_list = [50]
# elec_area = 5
# jpa_lns = 400
# jpa_lne = 450
# jpc_lns = 1310
# jpc_lne = 1330
# peak_range = 200
# peak_pos = 800
# trough_pos = 1600
# volt_arr1,current_arr1,scan_rate_sqrt1,jpa_arr1,jpc_arr1,lin_fit_jpa1,lin_fit_jpc1,peak_volt1,trough_volt1,jpa_base1,jpc_base1,peak_sep_list1 = cv_peak_sqrt(file_name1_1,file_name2_1,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.title('Ti', fontsize=16)
# plt.legend(fontsize=8)
# plt.xlabel('Potential (V)',fontsize=10)
# plt.ylabel('Current density (mA/cm$^2$)',fontsize=10)

# plt.show()
# plt.cla()
# plt.clf()


# fig, ax = plt.subplots(figsize=(15,15))
# file_name1_2 = 'CV and EIS/MX-UT0.1/mx-ut0.1-0to-1.2-'
# file_name2_2 = 'mvs.par'
# scan_list = [50,40,30,20,10,5,3]
# elec_area = 5
# jpa_lns = 400
# jpa_lne = 450
# jpc_lns = 1310
# jpc_lne = 1330
# peak_range = 200
# peak_pos = 800
# trough_pos = 1600
# volt_arr2,current_arr2,scan_rate_sqrt2,jpa_arr2,jpc_arr2,lin_fit_jpa2,lin_fit_jpc2,peak_volt2,trough_volt2,jpa_base2,jpc_base2,peak_sep_list2 = cv_peak_sqrt(file_name1_2,file_name2_2,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.title('MX-UT-CP-0.1mg/cm$^2$', fontsize=10)
# plt.legend(fontsize=8)
# plt.xlabel('Potential (V)',fontsize=10)
# plt.ylabel('Current density (mA/cm$^2$)',fontsize=10)
# plt.show()
# plt.cla()
# plt.clf()

# # fig, ax = plt.subplots(figsize=(15,15))
# file_name1_3 = 'CV and EIS/MX-HT-1mgcm/mc3_cv_MX_HT_1mgcm2_0to1.2_'
# file_name2_3 = 'mvs.par'
# scan_list = [50,40,30,20,10,5,3]
# elec_area = 5
# jpa_lns = 250
# jpa_lne = 450
# jpc_lns = 1310
# jpc_lne = 1330
# peak_range = 200
# peak_pos = 700
# trough_pos = 1600
# volt_arr3,current_arr3,scan_rate_sqrt3,jpa_arr3,jpc_arr3,lin_fit_jpa3,lin_fit_jpc3,peak_volt3,trough_volt3,jpa_base3,jpc_base3,peak_sep_list3 = cv_peak_sqrt(file_name1_3,file_name2_3,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.title('MX-HT-1mg/cm2', fontsize=16)
# plt.legend(fontsize=8)
# plt.xlabel('Potential (V)',fontsize=10)
# plt.ylabel('Current density (mA/cm$^2$)',fontsize=10)
# plt.show()
# plt.cla()
# plt.clf()

# fig, ax = plt.subplots(figsize=(15,15))
file_name1_4 = 'CV and EIS/MX-HT-1mgcm/mc3_cv_HT_1mgcm2_0to1.2_'
file_name2_4 = 'mvs.par'
scan_list = [50,40,30,20,10,5,3]
scan_list = [50]
elec_area = 5
jpa_lns = 300
jpa_lne = 450
jpc_lns = 1330
jpc_lne = 1350
peak_range = 200
peak_pos = 550
trough_pos = 1600
volt_arr4,current_arr4,scan_rate_sqrt4,jpa_arr4,jpc_arr4,lin_fit_jpa4,lin_fit_jpc4,peak_volt4,trough_volt4,jpa_base4,jpc_base4,peak_sep_list4 = cv_peak_sqrt(file_name1_4,file_name2_4,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
plt.title('HT', fontsize=10)
plt.legend(fontsize=8)
plt.xlabel('Potential (V)',fontsize=10)
plt.ylabel('Current density (mA/cm$^2$)',fontsize=10)
plt.show()
plt.cla()
plt.clf()

# # fig, ax = plt.subplots(figsize=(15,15))
# file_name1_5 = 'CV and EIS/wetted-ace-UT/UT-'
# file_name2_5 = 'mvs-0to1.2.par'
# scan_list = [50,30,20,10,5,3]
# elec_area = 5
# jpa_lns = 300
# jpa_lne = 450
# jpc_lns = 1280
# jpc_lne = 1330
# peak_range = 200
# peak_pos = 800
# trough_pos = 1600
# volt_arr5,current_arr5,scan_rate_sqrt5,jpa_arr5,jpc_arr5,lin_fit_jpa5,lin_fit_jpc5,peak_volt5,trough_volt5,jpa_base5,jpc_base4,peak_sep_list5 = cv_peak_sqrt(file_name1_5,file_name2_5,scan_list,elec_area,jpa_lns,jpa_lne,jpc_lns,jpc_lne,peak_range, peak_pos, trough_pos)
# plt.xlim(0,1.2)
# plt.title('UT-CP', fontsize=10)
# plt.legend(fontsize=8)
# plt.xlabel('Potential (V)',fontsize=10)
# plt.ylabel('Current density (mA/cm$^2$)',fontsize=10)
# plt.show()
# plt.cla()
# plt.clf()



# # fig, ax = plt.subplots(figsize=(10,10))
# scan_list = [50,40,30,20,10,5,3]
# plt.plot(np.sqrt(np.flip(scan_list)),jpa_arr1,'bo-',label='Ti')
# plt.plot(np.sqrt(np.flip(scan_list)),-jpc_arr1,'bo-')

# plt.plot(np.sqrt(np.flip(scan_list)),jpa_arr2,'s-',color='orange',label='MX-UT')
# plt.plot(np.sqrt(np.flip(scan_list)),-jpc_arr2,'s-',color='orange')

# jpa_arr3 = [2.11861,2.53683,3.21579,3.8,4.2,4.4,4.6]

# plt.plot(np.sqrt(np.flip(scan_list)),jpa_arr3,'g^-',label='MX-HT')
# plt.plot(np.sqrt(np.flip(scan_list)),-jpc_arr3,'g^-')

# plt.plot(np.sqrt(np.flip(scan_list)),jpa_arr4,'r*-',label='HT')
# plt.plot(np.sqrt(np.flip(scan_list)),-jpc_arr4,'r*-')

# scan_list = [50,30,20,10,5,3]
# plt.plot(np.sqrt(np.flip(scan_list)),jpa_arr5,'+-',color='purple',label='UT')
# plt.plot(np.sqrt(np.flip(scan_list)),-jpc_arr5,'+-',color='purple')

# plt.legend(fontsize=15)
# plt.xlabel('Scan rate$^{1/2}$ (mV/s)$^{1/2}$')
# plt.ylabel('Peak current (mA)')
# plt.show()
# plt.cla()
# plt.clf()

# # Reversibility
# # fig, ax = plt.subplots(figsize=(10,10))
# scan_list = [50,40,30,20,10,5,3]
# plt.plot(scan_list,jpa_arr1/jpc_arr1,'bo-',label='Ti')

# plt.plot(scan_list,jpa_arr2/jpc_arr2,'s-',color='orange',label='MX-UT')

# jpa_arr3 = np.array([2.11861,2.53683,3.21579,3.8,4.2,4.4,4.6])
# # print(jpa_arr3)
# # print(jpc_arr3.T)
# print(jpa_arr3/jpc_arr3.T)
# print(scan_list)
# # rev3 = jpa_arr3/jpc_arr3.T
# # print(rev3)
# plt.plot(scan_list,(jpa_arr3/jpc_arr3.T).T,'g^-',label='MX-HT')

# plt.plot(scan_list,jpa_arr4/jpc_arr4,'r*-',label='HT')


# scan_list = [50,30,20,10,5,3]
# plt.plot(scan_list,jpa_arr5/jpc_arr5,'+-',color='purple',label='UT')

# plt.legend(fontsize=15)
# plt.xlabel('Scan rate (mV/s)')
# plt.ylabel('jpa/jpc')
# plt.show()
# plt.cla()
# plt.clf()

# # fig, ax = plt.subplots(figsize=(10,10))
# scan_list = np.array([50,40,30,20,10,5,3])

# plt.plot(scan_list,peak_sep_list1*1000,'o-',label='Ti')
# plt.plot(scan_list[2:7],peak_sep_list2[2:7]*1000,'s-',label='MX-UT')
# plt.plot(scan_list,peak_sep_list3*1000,'^-',label='MX-HT')

# peak_sep_list4= np.array([0.663,0.611,0.532,0.463126,0.352099,0.267755,0.220829])
# plt.plot(scan_list,peak_sep_list4*1000,'*-',label='HT')

# scan_list = np.array([30,20,10,5,3])
# plt.plot(scan_list,peak_sep_list5[1:6]*1000,'+-',label='UT')

# # plt.plot(peak_sep_list3[1:3],scan_list[1:3],'^-',label='HT')
# # ht_x_3 = -0.4247881 - -0.6707667
# # ht_x_5 = -0.4165071 - -0.6913159
# # plt.plot([ht_x_5,ht_x_3],scan_list[1:3],'^-',label='HT')
# plt.legend(fontsize=15)
# plt.xlabel('Scan rate (mV/s)')
# plt.ylabel('Peak separation (mV)')
# # plt.xlim(0,1)
# # plt.ylim(0,50)
# plt.show()

# #EIS ##############################################

# fig, ax = plt.subplots(figsize=(3,3))

# eis_file3 = 'CV and EIS/MX-HT-1mgcm/eis_MX_UT_50vrms_0.67v.par'
# frequencies3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
# frequencies3, z3 = preprocessing.ignoreBelowX(frequencies3, z3)
# plot_nyquist(z3[:15],ax=ax,label='MX-UT')

# eis_file3 = 'CV and EIS/MX-HT-1mgcm/eis_HT_50vrms_0.67v.par'
# frequencies3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
# frequencies3, z3 = preprocessing.ignoreBelowX(frequencies3, z3)
# plot_nyquist(z3,ax=ax,label='HT')

# eis_file3 = 'CV and EIS/MX-HT-1mgcm/eis_MX_HT_50vrms_0.67v.par'
# frequencies3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
# frequencies3, z3 = preprocessing.ignoreBelowX(frequencies3, z3)
# plot_nyquist(z3,ax=ax,label='MX-HT')

# # eis_file3 = 'CV and EIS/EIS_new_met/eis_50mMv3_62.5mMh2so4_HT_50mVrms_-0.50Vref.par'
# # frequencies3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
# # frequencies3, z3 = preprocessing.ignoreBelowX(frequencies3, z3)
# # plot_nyquist(z3,ax=ax)

# plt.legend(fontsize=15)
# plt.show()
# plt.cla()
# plt.clf()

# #EIS ##############################################