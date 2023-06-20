import numpy as np 
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PySimpleCV_main_func import battery_xls2df, get_CV_init, find_state_seq, get_battery_eff, cy_idx_state_range, CV_file2df, get_CV_peak, eis_fit, eis_read_file, search_pattern


mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16 
fig, ax = plt.subplots(2, 2, figsize=(10,10))

color_lst = ["lightgreen","red","lightblue","orange","orchid"]
mark = ["o","s","^","v","D"]
bat_file = ["UT|HT","HT|HT", "MX-0.1|HT",  "MX-0.1|MX-0.1",  "MX-0.5|MX-0.5"]
col_idx = 0
for i in bat_file:
    file = "battery_ods/" + str(i) +".ods"
    df = pd.read_excel(file, engine="odf")
    ax[0,0].plot(df.loc[1:,'cycle number'],df.loc[1:,'charge-capacity'], label=str(i), linestyle='', marker=mark[col_idx], markerfacecolor=color_lst[col_idx], markeredgecolor='black',markeredgewidth=0.5,markersize=5)
    col_idx += 1
ax[0,0].legend(fontsize=11,loc = "upper left", bbox_to_anchor=(0.06,1))
ax[0,0].set_xlim(1,100)
ax[0,0].set_xlabel("Cycle Number")
ax[0,0].set_ylabel("Charge Capacity (mAh)")
ax[0,0].set_title('a)',y=1,x=0.05,pad=-14)    

col_idx = 0
for i in bat_file:
    file = "battery_ods/" + str(i) +".ods"
    df = pd.read_excel(file, engine="odf")
    ax[0,1].plot(df.loc[1:,'cycle number'],df.loc[1:,'discharge-capacity'], label=str(i), linestyle='', marker=mark[col_idx], markerfacecolor=color_lst[col_idx], markeredgecolor='black',markeredgewidth=0.5,markersize=5)
    dis_end = (df.loc[:,'discharge-capacity']).to_numpy()[-1]
    dis_max = max((df.loc[:,'discharge-capacity']).to_numpy())
    dis_start = ((df.loc[:,'discharge-capacity']).to_numpy())[0]
    print(str(i),"DC At 100 cycle",dis_end)
    # print(str(i),"Max DC Cap",dis_max)
    print(str(i),"DC cap ret",dis_end*100/dis_max)
    print(str(i),"DC cap ret start",dis_end*100/dis_start)
    col_idx += 1
# ax[0,1].legend(fontsize=8)
ax[0,1].set_xlim(1,100)
ax[0,1].set_xlabel("Cycle Number")
ax[0,1].set_ylabel("Discharge Capacity (mAh)")
ax[0,1].set_title('b)',y=1,x=0.95,pad=-14)    

col_idx = 0
for i in bat_file:
    file = "battery_ods/" + str(i) +".ods"
    df = pd.read_excel(file, engine="odf")
    ax[1,0].plot(df.loc[1:,'cycle number'],df.loc[1:,'ee'], label=str(i), linestyle='', marker=mark[col_idx], markerfacecolor=color_lst[col_idx], markeredgecolor='black',markeredgewidth=0.5,markersize=5)
    col_idx += 1
# ax[1,0].legend(fontsize=8)
ax[1,0].set_xlim(1,100)
ax[1,0].set_ylim(0,100)
ax[1,0].set_xlabel("Cycle Number")
ax[1,0].set_ylabel("Energy Efficiency (%)")
ax[1,0].set_title('c)',y=1,x=0.05,pad=-14)    

categories = ["Voltage", "Coulombic", "Energy"]
# Width of each bar
bar_width = 0.2

# Generate an array for x-values for each category
x = np.arange(len(categories))

# Shift the bars
x_bar1 = x - bar_width
x_bar2 = x
x_bar3 = x + bar_width

bat_file = ["UT|HT","HT|HT", "MX-0.1|HT",  "MX-0.1|MX-0.1",  "MX-0.5|MX-0.5"]
utht = []
wid = -0.3
col_idx = 0
for i in bat_file:
    file = "battery_ods/" + str(i) +".ods"
    df = pd.read_excel(file, engine="odf")
    avg_ve = np.average(df.loc[1:,'ve'])
    avg_ce = np.average(df.loc[1:,'ce'])
    avg_ee = np.average(df.loc[1:,'ee'])   
    elec = [avg_ve,avg_ce,avg_ee]   
    rect = ax[1,1].bar(x+wid, elec, width=0.15, label=str(i), color=color_lst[col_idx])
    ax[1,1].bar_label(rect, padding=-12, fmt='%.0f',color='black')
    wid += 0.15
    col_idx += 1

categories = ["Voltage", "Coulombic",  "Energy"]
ax[1,1].set_ylabel("Efficiencies (%)")
ax[1,1].set_xticks(x, categories)
ax[1,1].set_title('d)',y=1,x=0.05,pad=-14)   
ax[1,1].set_ylim(0,100)
# Display the graph
plt.show()
















# fig, ax = plt.subplots(figsize=(13,13))
# mx_ut_cap = df.iloc[2:, 1:3]
# # plt.plot(mx_ut_cap.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='lightblue', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# plt.plot(mx_ut_cap.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='lightblue', markeredgecolor='black',markeredgewidth=0.5,markersize=5)

# ht_cap = df.iloc[2:, 8:10]
# # plt.plot(ht_cap.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='red', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# plt.plot(ht_cap.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='red', markeredgecolor='black',markeredgewidth=0.5,markersize=5)

# ut_cap = df.iloc[2:, 15:17]
# # plt.plot(ut_cap.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='lightgreen', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# plt.plot(ut_cap.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='lightgreen', markeredgecolor='black',markeredgewidth=0.5,markersize=5)

# mxmx_cap = df.iloc[2:, 22:24]
# # plt.plot(mxmx_cap.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='orange', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# plt.plot(mxmx_cap.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='orange', markeredgecolor='black',markeredgewidth=0.5,markersize=5)

# plt.xlabel("Cycle number")
# plt.ylabel("Capacity (mAh)")
# plt.xlim(0,100)
# plt.ylim(0,200)

# # cap_leg = Line2D([0], [0], label='Charge capacity', markerfacecolor='None', markeredgecolor='black',markeredgewidth=1, marker='^', linestyle='None')
# # discap_leg = Line2D([0], [0], label='Discharge capacity', markerfacecolor='None', markeredgecolor='black',markeredgewidth=1, marker='o', linestyle='None')
# lab_leg = Patch(facecolor='None', label='Anode | Cathode')
# mx_ut_leg = Patch(facecolor='lightblue', label='MX-UT | HT')
# ht_leg = Patch(facecolor='red', label='HT | HT')
# ut_leg = Patch(facecolor='lightgreen', label='UT | HT')
# mxmx_leg = Patch(facecolor='orange', label='MX-UT | MX-UT')

# plt.legend(handles=[lab_leg,mx_ut_leg,mxmx_leg,ht_leg,ut_leg], loc='upper right',fontsize="7")
# plt.show()

# mx_ut_eff = df.iloc[2:, 3:6]
# # plt.plot(mx_ut_eff.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='lightblue', markeredgecolor='black',markeredgewidth=0.5)
# # plt.plot(mx_ut_eff.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='lightblue', markeredgecolor='black',markeredgewidth=0.5)
# plt.plot(mx_ut_eff.iloc[:,2:3], linestyle='', marker='o', markerfacecolor='lightblue', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# avg_mx_ut_VE = np.average(mx_ut_eff.iloc[:,0:1])
# avg_mx_ut_CE = np.average(mx_ut_eff.iloc[:,1:2])
# avg_mx_ut_EE = np.average(mx_ut_eff.iloc[:,2:3])

# ht_eff = df.iloc[2:, 10:13]
# # plt.plot(ht_eff.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='red', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# # plt.plot(ht_eff.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='red', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# plt.plot(ht_eff.iloc[:,2:3], linestyle='', marker='o', markerfacecolor='red', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# avg_ht_VE = np.average(ht_eff.iloc[:,0:1])
# avg_ht_CE = np.average(ht_eff.iloc[:,1:2])
# avg_ht_EE = np.average(ht_eff.iloc[:,2:3])

# ut_eff = df.iloc[2:, 17:20]
# # plt.plot(ut_eff.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='lightgreen', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# # plt.plot(ut_eff.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='lightgreen', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# plt.plot(ut_eff.iloc[:,2:3], linestyle='', marker='o', markerfacecolor='lightgreen', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# avg_ut_VE = np.average(ut_eff.iloc[:,0:1])
# avg_ut_CE = np.average(ut_eff.iloc[:,1:2])
# avg_ut_EE = np.average(ut_eff.iloc[:,2:3])
# # print(ut_eff)

# mxmx_eff = df.iloc[2:, 24:27]
# # plt.plot(mxmx_eff.iloc[:,0:1], linestyle='', marker='^', markerfacecolor='orange', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# # plt.plot(mxmx_eff.iloc[:,1:2], linestyle='', marker='o', markerfacecolor='orange', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# plt.plot(mxmx_eff.iloc[:,2:3], linestyle='', marker='o', markerfacecolor='orange', markeredgecolor='black',markeredgewidth=0.5,markersize=5)
# avg_mxmx_VE = np.average(mxmx_eff.iloc[:,0:1])
# avg_mxmx_CE = np.average(mxmx_eff.iloc[:,1:2])
# avg_mxmx_EE = np.average(mxmx_eff.iloc[:,2:3])

# plt.xlabel("Cycle number")
# plt.ylabel("Efficiency (%)")
# plt.xlim(0,100)
# plt.ylim(0,100)

# # VE_leg = Line2D([0], [0], label='Voltage Efficiency', markerfacecolor='None', markeredgecolor='black',markeredgewidth=1, marker='^', linestyle='None')
# # CE_leg = Line2D([0], [0], label='Coulombic Efficiency', markerfacecolor='None', markeredgecolor='black',markeredgewidth=1, marker='o', linestyle='None')
# EE_leg = Line2D([0], [0], label='Energy Efficiency', markerfacecolor='None', markeredgecolor='black',markeredgewidth=1, marker='s', linestyle='None')
# lab_leg = Patch(facecolor='None', label='Anode | Cathode')
# mx_ut_leg = Patch(facecolor='lightblue', label='MX-UT | HT')
# mxmx_leg = Patch(facecolor='orange', label='MX-UT | MX-UT')
# ht_leg = Patch(facecolor='red', label='HT | HT')
# ut_leg = Patch(facecolor='lightgreen', label='UT | HT')

# plt.legend(handles=[lab_leg,mx_ut_leg,mxmx_leg,ht_leg,ut_leg],loc='lower left',fontsize="7")
# plt.show()

# # print(avg_ut_VE)
# species = ("Voltage", "Coulombic", "Energy")
# electrodes = {
#     'MX-UT | HT': (round(avg_mx_ut_VE,0), round(avg_mx_ut_CE,0), round(avg_mx_ut_EE,0)),
#     'MX-UT | MX-UT': (round(avg_mxmx_VE,0), round(avg_mxmx_CE,0), round(avg_mxmx_EE,0)),
#     'HT | HT': (round(avg_ht_VE,0), round(avg_ht_CE,0), round(avg_ht_EE,0)),
#     'UT | HT': (round(avg_ut_VE,0), round(avg_ut_CE,0), round(avg_ut_EE,0)),
# }

# x = np.arange(len(species))  # the label locations
# width = 0.2  # the width of the bars
# multiplier = 0

# fig, ax = plt.subplots(figsize=(3,4),layout='constrained')
# bar_color = ['lightblue','orange', 'red', 'lightgreen']
# for attribute, measurement in electrodes.items():
#     offset = width * multiplier
#     rects = ax.bar(x + offset, measurement, width, label=attribute, color=bar_color[multiplier])
#     ax.bar_label(rects, padding=-12)
#     multiplier += 1

# ax.set_ylabel('Efficiency (%)')
# ax.set_xticks(x + width, species)
# # handles = ['Anode | Cathode','MX-UT-CP 0.1mg/cm$^2$ | HT-CP', 'HT-CP | HT-CP', 'UT-CP | HT-CP']
# # ax.legend(handles=[lab_leg,mx_ut_leg,ht_leg,ut_leg],loc='upper left', ncols=1,fontsize="9")
# ax.set_ylim(0, 100)

# plt.show()



