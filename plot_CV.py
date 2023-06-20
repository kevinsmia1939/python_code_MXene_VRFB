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

cv_file = [
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_50mVs_beaker.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_40mVs_beaker.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_30mVs_beaker.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_20mVs_beaker.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_10mVs_beaker_4.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_5mVs_beaker.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_3mVs_beaker.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_3mVs_beaker_4.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_3mVs_beaker_3.par',
            
            'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_1mVs_beaker_3.par',
            'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_1mVs_beaker_4.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_50mVs_beaker.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to1.2V_30mVs_beaker_2.par',
            # 'cv_v4_0.05MH2SO4_50mMV4_HTSig_0to-0.8V_40mVs_beaker.par'
            ]



plt.figure(figsize=(10,10))

cut_val = 20
for i in cv_file:
    df_CV = CV_file2df('CV and EIS/low_h2so4_conc/'+i)
    cv_size, volt, current = get_CV_init(df_CV)
    plt.plot(volt[cut_val:], current[cut_val:]/5)

plt.xlabel('Voltage (V) vs Ag/AgCl', fontsize=16)
plt.ylabel('Current density (A/cm$^2$)', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlim(-1, 1.5)
# plt.ylim(-0.005, 0.003)
# plt.grid()
plt.show()

# fig_eis = Figure(figsize=(10,10))
# ax_eis = fig_eis.add_subplot(111)
# ax_eis.set_xlabel("Z'(w)[Ohms]")
# ax_eis.set_ylabel('-Z"(w)[Ohms]')
# ax_eis.set_aspect('equal', 'box')
# ax_eis.axis('equal')
# ax_eis.grid()

# # plt.figure(figsize=(10,10))
# fig, ax = plt.subplots(figsize=(10,10))
# eis_file = ['eis_v4_0.5M_HT_0.8vvsref_sigracet_5cm2.par',
#             'eis_v4_0.5M_HT_0.84vvsref_sigracet_5cm2.par',
#             'eis_v4_0.5M_HT_0.92vvsref_sigracet_5cm2.par',
#             'eis_v4_0.5M_HT_0.96vvsref_sigracet_5cm2.par',
#             'eis_v4_0.5M_HT_0.884vvsref_sigracet_5cm2.par',
#             ]
# for i in eis_file:
#     eis_file = 'CV and EIS/'+i
#     frequencies, z = preprocessing.readFile(eis_file,instrument='versastudio')
#     frequencies, z = preprocessing.ignoreBelowX(frequencies, z) 
#     plot_nyquist(ax, z,fmt='-')
#     # circuit.plot(f_data=frequencies, Z_data=z, kind='nyquist')


# plt.show()

# # from impedance import preprocessing
# # from impedance.models.circuits import CustomCircuit
# # import matplotlib.pyplot as plt
# # from impedance.visualization import plot_nyquist
# # import numpy as np
# # from impedance.models.circuits import Randles, CustomCircuit
# # eis_file = 'CV and EIS/eis_v4_0.5M_HT_0.884vvsref_sigracet_5cm2.par'
# # frequencies, z = preprocessing.readVersaStudio(eis_file)
# # frequencies, z = preprocessing.ignoreBelowX(frequencies, z)

# # fig, ax = plt.subplots(figsize=(10,10))
# # # circuit = CustomCircuit('R0-p(R1-W1,CPE1)', initial_guess=[1.2],constants={'R1':0.35, 'CPE1_0':0.7, 'CPE1_1':0.6,'W1':0.07})
# # # circuit = CustomCircuit('R0-p(R1-W1,CPE1)', initial_guess=[1.2,0.35,0.7,0.6,0.07])
# # # circuit = Randles(initial_guess=[1.2,0.35,0.7,0.6,0.07,0.5],CPE=True)
# # # f_pred = np.logspace(5,2)
# # # circuit.fit(frequencies, z)
# # # circuit_fit = circuit.predict(f_pred)
# # plot_nyquist(ax, z,fmt='-')
# # # plot_nyquist(ax, circuit_fit, fmt='--')
# # plt.show()
# # # print(circuit)