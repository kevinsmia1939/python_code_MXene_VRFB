import numpy as np 
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredText
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

cut_val = 0
elec_area = 10
fig, ax = plt.subplots(3, 2, figsize=(20,20))


df_CV1 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_1mVs_beaker.par')
_, volt1, current1 = get_CV_init(df_CV1)
ax[0, 0].plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='1 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_3mVs_beaker.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[0, 0].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='3 mV/s')
df_CV3 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_5mVs_beaker.par')
_, volt3, current3 = get_CV_init(df_CV3)
ax[0, 0].plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='5 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_10mVs_beaker.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[0, 0].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='10 mV/s')

# v_lin = np.linspace(0,-0.8,30)
# print(v_lin)
# ax[0, 0].plot(v_lin, (v_lin*0.2)**3,label='10 mV/s')

ax[0, 0].set_xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
ax[0, 0].set_ylabel('Current density (A/cm$^2$)', fontsize=16)
ax[0, 0].legend(fontsize=15)
ax[0, 0].set_xlim(-0.8, 0)
ax[0, 0].set_ylim(-0.005, 0.002)


df_CV1 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to1.2V_1mVs_beaker.par')
_, volt1, current1 = get_CV_init(df_CV1)
ax[0, 1].plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='1 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to1.2V_3mVs_beaker.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[0, 1].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='3 mV/s')
df_CV3 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to1.2V_5mVs_beaker.par')
_, volt3, current3 = get_CV_init(df_CV3)
ax[0, 1].plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='5 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to1.2V_10mVs_beaker.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[0, 1].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='10 mV/s')
ax[0, 1].set_xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
ax[0, 1].set_ylabel('Current density (A/cm$^2$)', fontsize=16)
ax[0, 1].legend(fontsize=15)
ax[0, 1].set_xlim(0, 1.2)
ax[0, 1].set_ylim(-0.004, 0.004)

df_CV1 = CV_file2df('CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_1mVs_beaker_2.par')
_, volt1, current1 = get_CV_init(df_CV1)
ax[1, 1].plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='1 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_3mVs_beaker_2.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[1, 1].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='3 mV/s')
df_CV3 = CV_file2df('CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_5mVs_beaker_2.par')
_, volt3, current3 = get_CV_init(df_CV3)
ax[1, 1].plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='5 mV/s')
df_CV4 = CV_file2df('CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_10mVs_beaker_2.par')
_, volt4, current4 = get_CV_init(df_CV4)
ax[1, 1].plot(volt4[cut_val:], current4[cut_val:]/elec_area,label='10 mV/s')

df_CV5 = CV_file2df('CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_20mVs_beaker_2.par')
_, volt5, current5 = get_CV_init(df_CV5)
ax[1, 1].plot(volt5[cut_val:], current5[cut_val:]/elec_area,label='20 mV/s')
df_CV6 = CV_file2df('CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_30mVs_beaker_2.par')
_, volt6, current6 = get_CV_init(df_CV6)
ax[1, 1].plot(volt6[cut_val:], current6[cut_val:]/elec_area,label='30 mV/s')

df_CV7 = CV_file2df('CV and EIS/low_h2so4_conc/cv_v4v5_HTsig_2/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_50mVs_beaker_2.par')
_, volt7, current7 = get_CV_init(df_CV7)
ax[1, 1].plot(volt7[cut_val:], current7[cut_val:]/elec_area,label='50 mV/s')

ax[1, 1].set_xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
ax[1, 1].set_ylabel('Current density (A/cm$^2$)', fontsize=16)
ax[1, 1].legend(fontsize=15)
ax[1, 1].set_xlim(0, 1.2)
# ax[1, 1].set_ylim(-0.004, 0.0044)

df_CV1 = CV_file2df('CV and EIS/low_h2so4_conc/sigracet_v2v3/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_1mVs_beaker.par')
_, volt1, current1 = get_CV_init(df_CV1)
ax[1, 0].plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='1 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/sigracet_v2v3/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_3mVs_beaker_4.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[1, 0].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='3 mV/s')
df_CV3 = CV_file2df('CV and EIS/low_h2so4_conc/sigracet_v2v3/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_5mVs_beaker.par')
_, volt3, current3 = get_CV_init(df_CV3)
ax[1, 0].plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='5 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/sigracet_v2v3/cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_10mVs_beaker_4.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[1, 0].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='10 mV/s')
ax[1, 0].set_xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
ax[1, 0].set_ylabel('Current density (A/cm$^2$)', fontsize=16)
ax[1, 0].legend(fontsize=15)
ax[1, 0].set_xlim(-0.8,0)
ax[1, 0].set_ylim(-0.005, 0.002)

df_CV1 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to-0.65v_1mvs_0.1mgcm-mx_test.par')
_, volt1, current1 = get_CV_init(df_CV1)
ax[2, 0].plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='1 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to-0.8v_3mvs_0.3mgcm-mx_test.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[2, 0].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='3 mV/s')
df_CV3 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to-0.8v_5mvs_0.3mgcm-mx_test.par')
_, volt3, current3 = get_CV_init(df_CV3)
ax[2, 0].plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='5 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to-0.8v_10mvs_0.3mgcm-mx_test.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[2, 0].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='10 mV/s')
ax[2, 0].set_xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
ax[2, 0].set_ylabel('Current density (A/cm$^2$)', fontsize=16)
ax[2, 0].legend(fontsize=15)
ax[2, 0].set_xlim(-0.8,0)
ax[2, 0].set_ylim(-0.005, 0.002)

df_CV1 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to1.2v_1mvs_0.1mgcm-mx_test.par')
_, volt1, current1 = get_CV_init(df_CV1)
ax[2, 1].plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='1 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to1.2v_3mvs_0.3mgcm-mx_test.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[2, 1].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='3 mV/s')
df_CV3 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to1.2v_5mvs_0.3mgcm-mx_test.par')
_, volt3, current3 = get_CV_init(df_CV3)
ax[2, 1].plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='5 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/mxene/mx_v4_50mMh2so4_50mMV4_HTSig_0to1.2v_10mvs_0.3mgcm-mx_test.par')
_, volt2, current2 = get_CV_init(df_CV2)
ax[2, 1].plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='10 mV/s')
ax[2, 1].set_xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
ax[2, 1].set_ylabel('Current density (A/cm$^2$)', fontsize=16)
ax[2, 1].legend(fontsize=15)
ax[2, 1].set_xlim(0,1.2)
ax[2, 1].set_ylim(-0.004, 0.004)

at1 = AnchoredText("Pristine carbon paper", prop=dict(size=15), frameon=True, loc='upper left')
at1.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax[0,0].add_artist(at1)

at2 = AnchoredText("Pristine carbon paper", prop=dict(size=15), frameon=True, loc='upper right')
at2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax[0,1].add_artist(at2)

at3 = AnchoredText("Heat-treated carbon paper", prop=dict(size=15), frameon=True, loc='upper right')
at3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax[1,1].add_artist(at3)

at4 = AnchoredText("Heat-treated carbon paper", prop=dict(size=13), frameon=True, loc='upper right')
at4.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax[1,0].add_artist(at4)

at5 = AnchoredText("MXene 0.3 mg/cm2", prop=dict(size=15), frameon=True, loc='upper right')
at5.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax[2,0].add_artist(at5)

at6 = AnchoredText("MXene 0.3 mg/cm2", prop=dict(size=15), frameon=True, loc='upper right')
at6.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax[2,1].add_artist(at6)
# loc = plticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
# ax.xaxis.set_major_locator(loc)
plt.show()
