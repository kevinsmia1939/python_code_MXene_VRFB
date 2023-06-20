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
fig, ax = plt.subplots(figsize=(15,15))


df_CV1 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_1mVs_beaker.par')
_, volt1, current1 = get_CV_init(df_CV1)
plt.plot(volt1[cut_val:], current1[cut_val:]/elec_area,label='1 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_3mVs_beaker.par')
_, volt2, current2 = get_CV_init(df_CV2)
plt.plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='3 mV/s')
df_CV3 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_5mVs_beaker.par')
_, volt3, current3 = get_CV_init(df_CV3)
plt.plot(volt3[cut_val:], current3[cut_val:]/elec_area,label='5 mV/s')
df_CV2 = CV_file2df('CV and EIS/low_h2so4_conc/untreated/cv_v4_0.05MH2SO4_50mMV4_UTSig_0to-0.8V_10mVs_beaker.par')
_, volt2, current2 = get_CV_init(df_CV2)
plt.plot(volt2[cut_val:], current2[cut_val:]/elec_area,label='10 mV/s')
plt.xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
plt.ylabel('Current density (A/cm$^2$)', fontsize=16)
plt.legend(fontsize=15)
plt.xlim(-0.8, 0)
plt.ylim(-0.005, 0.002)