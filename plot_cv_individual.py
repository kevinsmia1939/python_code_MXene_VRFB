import numpy as np 
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from impedance import preprocessing
from impedance.models.circuits import Randles, CustomCircuit
from impedance.visualization import plot_nyquist

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

cut_val = 3998


# fig, ax = plt.subplots(figsize=(20,10))
# eis_file3 = 'CV and EIS/low_h2so4_conc/mxene/EIS/eis_2Mh2so4_50mMV4_1mgcm2_10mV_0.7V.par'
# freq3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
# freq3, z3 = preprocessing.ignoreBelowX(freq3, z3)
# # circuit.plot(f_data=freq, Z_data=Z, kind='nyquist')
# plot_nyquist(z3,ax=ax)
# # ax.plot(np.real(z_fit), -np.imag(z_fit),'--')
# # ax.plot(np.real(z3), -np.imag(z3),':')
# # ax.legend(fontsize=15)
# plt.show()

# plt.cla()
# plt.clf()

# fig, ax = plt.subplots(figsize=(10,10))
# start_val = 3995
# end_val = 6000

# df_CV4 = CV_file2df('CV and EIS/wetted-ace-UT/UT-5mvs-0to-0.75.par')
# _, volt4, current4 = get_CV_init(df_CV4)
# plt.plot(volt4, current4/5*1000,label='Untreated CP')

# df_CV3 = CV_file2df('CV and EIS/cv_new_method_2/cv_mc_50mMv4_50mMh2so4_HT_0to-0.75_5mvs.par')
# _, volt3, current3 = get_CV_init(df_CV3)
# plt.plot(volt3[start_val:end_val], current3[start_val:end_val]/5*1000,label='Heat treated CP')

# df_CV1 = CV_file2df('CV and EIS/MX-UT_5mvs.par')
# _, volt1, current1 = get_CV_init(df_CV1)
# plt.plot(volt1[0:2000], current1[0:2000]/5*1000,label='MXene-Untreated CP 1mg/cm$^2$')

# df_CV2 = CV_file2df('CV and EIS/cv_new_method_2/cv_mc_50mMv4_50mMh2so4_MX-HT_0to-0.75_5mvs.par')
# _, volt2, current2 = get_CV_init(df_CV2)
# plt.plot(volt2[start_val:end_val], current2[start_val:end_val]/5*1000,label='MXene-Heat treated CP 1mg/cm$^2$')

# plt.xlim(-0.8,0)
# plt.xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
# plt.ylabel('Current density (mA/cm$^2$)', fontsize=16)
# plt.legend(fontsize=15)
# plt.ylim(-6,2)
# plt.show()
# plt.cla()
# plt.clf()

# fig, ax = plt.subplots(figsize=(10,10))
# start_val = 3995
# end_val = 6000



# df_CV2 = CV_file2df('CV and EIS/cv_new_method_2/cv_mc_50mMv4_50mMh2so4_MX-HT_0to1.2_5mvs.par')
# _, volt2, current2 = get_CV_init(df_CV2)
# plt.plot(volt2[start_val:end_val], current2[start_val:end_val]/5*1000,label='MX-HT 1mg/cm2 5 mV/s')

# df_CV3 = CV_file2df('CV and EIS/cv_new_method_2/cv_mc_50mMv4_50mMh2so4_HT_0to1.2_5mvs.par')
# _, volt3, current3 = get_CV_init(df_CV3)
# plt.plot(volt3[start_val:end_val], current3[start_val:end_val]/5*1000,label='HT 5 mV/s')

# plt.xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
# plt.ylabel('Current density (mA/cm$^2$)', fontsize=16)
# plt.legend(fontsize=15)
# plt.show()
# plt.cla()
# plt.clf()

#EIS ##############################################

fig, ax = plt.subplots(figsize=(15,10))

eis_file3 = 'CV and EIS/EIS_new_met/eis_50mMv3_62.5mMh2so4_HT_50mVrms_-0.53Vref.par'
freq3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
freq3, z3 = preprocessing.ignoreBelowX(freq3, z3)
freq3 = freq3[0:18]
z3 = z3[0:18]
circuit = CustomCircuit(initial_guess=[1.45, .27, 2.59e-03, 7.76e-01], circuit='R_0-p(R_1,CPE_1)')
circuit.fit(freq3, z3, weight_by_modulus=True, global_opt=False)
ax = circuit.plot(ax, f_data=freq3, Z_data=z3, kind='nyquist',label='HT')
print(circuit)

eis_file3 = 'CV and EIS/EIS_new_met/eis_50mMv3_62.5mMh2so4_MX0.3_50mVrms_-0.53Vref.par'
freq3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
freq3, z3 = preprocessing.ignoreBelowX(freq3, z3)
freq3 = freq3[0:18]
z3 = z3[0:18]
circuit = CustomCircuit(initial_guess=[1.45, .27, 2.59e-03, 7.76e-01], circuit='R_0-p(R_1,CPE_1)')
circuit.fit(freq3, z3, weight_by_modulus=True, global_opt=False)
ax = circuit.plot(ax, f_data=freq3, Z_data=z3, kind='nyquist',label='MX')
# ax.set_xlim(0)


# circuit.fit(freq, z, sigma=sigma)

# z_fit = circuit.predict(freq)

# fig, ax = plt.subplots(figsize=(10,10))
# ax = circuit.plot(ax, f_data=freq, Z_data=z, kind='nyquist')

# eis_file3 = 'CV and EIS/EIS_new_met/eis_50mMv3_62.5mMh2so4_HT_50mVrms_-0.51Vref.par'
# freq3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
# freq3, z3 = preprocessing.ignoreBelowX(freq3, z3)
# plot_nyquist(z3,ax=ax)

# eis_file3 = 'CV and EIS/EIS_new_met/eis_50mMv3_62.5mMh2so4_HT_50mVrms_-0.50Vref.par'
# freq3, z3 = preprocessing.readFile(eis_file3,instrument='versastudio')
# freq3, z3 = preprocessing.ignoreBelowX(freq3, z3)
# plot_nyquist(z3,ax=ax)

plt.legend(fontsize=15)
plt.show()
plt.cla()
plt.clf()

#EIS ##############################################

# fig, ax = plt.subplots(figsize=(10,10))
# start_val = 1995
# end_val = 4000

# df_CV1 = CV_file2df('CV and EIS/wetted-ace-UT/UT-5mvs-0to1.2.par')
# _, volt1, current1 = get_CV_init(df_CV1)
# plt.plot(volt1, current1/5*1000,label='Untreated CP')

# df_CV2 = CV_file2df('CV and EIS/cv_new_method_2/cv_mc_50mMv4_50mMh2so4_Ti_0to1.2_5mvs.par')
# _, volt2, current2 = get_CV_init(df_CV2)
# plt.plot(volt2[start_val:end_val], current2[start_val:end_val]/1.2*1000,label='Ti 5 mV/s')

# df_CV3 = CV_file2df('CV and EIS/cv_new_method_2/cv_mc_50mMv4_50mMh2so4_HT_0to1.2_5mvs.par')
# _, volt3, current3 = get_CV_init(df_CV3)
# plt.plot(volt3[start_val:end_val], current3[start_val:end_val]/5*1000,label='Heat treated')



# df_CV4 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to1.2v_5mvs_1mgcm-mx_test.par')
# _, volt4, current4 = get_CV_init(df_CV4)
# plt.plot(volt4, current4/5*1000,label='MXene-Heat treated CP 1mg/cm$^2$')

# curr_min = current3[start_val:end_val]
# aab = volt3[np.where(curr_min-min(curr_min)==0)[0]]
# # plt.plot(aab,volt3[])
# plt.vlines(aab,-1,-4,linestyle='--')

# plt.xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
# plt.ylabel('Current density (mA/cm$^2$)', fontsize=16)
# plt.legend(fontsize=15)
# plt.xlim(0,1.2)
# plt.ylim(-4,7)
# plt.show()
# plt.cla()
# plt.clf()

fig, ax = plt.subplots(figsize=(10,10))
start_val = 1995
end_val = 4000
df_CV2 = CV_file2df('CV and EIS/cv_new_method_2/cv_mc_50mMv4_50mMh2so4_Ti_0to-0.75_5mvs.par')
_, volt2, current2 = get_CV_init(df_CV2)
plt.plot(volt2[start_val:end_val], current2[start_val:end_val]/1.2*1000,label='Ti 5 mV/s')

plt.xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
plt.ylabel('Current density (mA/cm$^2$)', fontsize=16)
plt.legend(fontsize=15)
plt.show()
plt.cla()
plt.clf()

fig, ax = plt.subplots(figsize=(10,10))
elec_area = 5*1000
df_CV1 = CV_file2df('CV and EIS/MX-HT-1mgcm/mc3_cv_Ti_0to1.2_5mvs.par')
_, volt1, current1 = get_CV_init(df_CV1)
plt.plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='Ti 5 mV/s')
df_CV2 = CV_file2df('CV and EIS/MX-HT-1mgcm/mc3_cv_MX_UT_1mgcm2_0to1.2_5mvs.par')
_, volt2, current2 = get_CV_init(df_CV2)
plt.plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='MX-HT')
df_CV3 = CV_file2df('CV and EIS/MX-HT-1mgcm/mc3_cv_MX_HT_1mgcm2_0to1.2_5mvs.par')
_, volt3, current3 = get_CV_init(df_CV3)
plt.plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='HT')
df_CV4 = CV_file2df('CV and EIS/MX-HT-1mgcm/mc3_cv_HT_1mgcm2_0to1.2_5mvs.par')
_, volt4, current4 = get_CV_init(df_CV4)
plt.plot(volt4[cut_val:], current4[cut_val:]/elec_area,label='10 mV/s')

df_CV5 = CV_file2df('CV and EIS/wetted-ace-UT/UT-5mvs-0to1.2.par')
_, volt5, current5 = get_CV_init(df_CV5)
plt.plot(volt5[cut_val:], current5[cut_val:]/elec_area,label='UT')
plt.xlim(0,1.2)
plt.legend()
plt.show()
plt.cla()
plt.clf()