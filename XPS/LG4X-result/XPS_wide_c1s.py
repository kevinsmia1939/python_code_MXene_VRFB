import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 70
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.titlesize'] = 25
mpl.rcParams['axes.labelsize'] = 16 
fig, ax = plt.subplots(figsize=(20,15))

label_size = 24

ax1 = plt.subplot(121)
ax1.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax2 = plt.subplot(122)
ax2.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax2.set_xlabel('Binding energy (eV)')
ax2.set_ylabel('Intensity (arb. unit)')


xps_lst = ['UT','HT','MX-UT-0.1','MX-UT-0.5']
xps_txt = ['Untreated carbon paper','Heat-treated carbon paper','MXene-coated \n (0.1mg/cm$^2$)','MXene-coated \n (0.5mg/cm$^2$)']
color_lst = ['green','red','royalblue','orange']
shift_idx = 0
shift=[-44000,-40000,-20000,0]
shift_l = 0

ax1.vlines(533,-350000,100000,colors='black',linestyles='--')
ax1.vlines(685.2,-350000,-100000,colors='black',linestyles='--')
ax1.vlines(200,-350000,-190000,colors='black',linestyles='--')
ax1.vlines(456,-350000,-100000,colors='black',linestyles='--')

for xps in xps_lst:
    # print(shift[shift_idx])
    excel = '../Excel sheets/' + xps + '.xlsx'
    df = pd.read_excel(excel,sheet_name='0-wide')  
    bind = df.iloc[1:,0].to_numpy() #Put C1s to 284.5 eV
    count = df.iloc[1:,1].to_numpy()-shift_l#+shift[shift_idx]
    ax1.plot(bind,count,linewidth = 2,color=color_lst[shift_idx])
    ax1.text(950, -shift_l+140000+shift[shift_idx], xps_txt[shift_idx], rotation=0, fontsize=25,horizontalalignment='center',verticalalignment='center')
    shift_l += 100000
    shift_idx += 1

ax1.set_xlabel('Binding energy (eV)')
ax1.set_ylabel('Intensity (arb. unit)')
ax1.set_xlim(bind[0],bind[-1])
ax1.set_ylim(-305000,440000)
ax1.text(284.5, 431000, 'C 1s', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
ax1.text(450, -60000, 'Ti 2p', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
# ax1.text(450, -139000, 'Ti 2p', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
ax1.text(525, 112000, 'O 1s', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
# ax1.text(525, 11000, 'O 1s', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
# ax1.text(525, -80000, 'O 1s', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
# ax1.text(525, -165000, 'O 1s', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
ax1.text(200, -175000, 'Cl 2p', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
# ax1.text(200, -275000, 'Cl 2p', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
ax1.text(685, -85000, 'F 1s', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')
# ax1.text(685, -165000, 'F 1s', rotation=0, fontsize=label_size,horizontalalignment='center',verticalalignment='center')

ax1.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)

def plot_xps(file_path,axis,mv_res,shift):
    df = pd.read_csv(file_path+'.csv', delimiter=',', header=1)
    
    bind = df['corrected x'].to_numpy()
    count = df['raw_y'].to_numpy()
    sum_fit = df['sum_fit'].to_numpy()
    residual = count - sum_fit
    std_res = np.std(residual)
    bg = df['bg'].to_numpy()
    print("##############################################")
    with open(file_path+'.txt', 'r') as file:
        lines = file.readlines()
        peaks = lines[11].split()
    print(peaks)
    # print(df.columns.tolist())
    df = df.drop('raw_x', axis=1)
    df = df.drop('raw_y', axis=1)
    df = df.drop('corrected x', axis=1)
    df = df.drop('data-bg', axis=1)
    df = df.drop('sum_components', axis=1)
    df = df.drop('bg', axis=1)
    df = df.drop('sum_fit', axis=1)
    try:
        df = df.drop('bg_shirley_', axis=1)
    except:
        pass
    size = df.shape[1]
    for i in np.arange(0,size):
        header = str(df.columns.tolist()[i])
        # print(header)
        # print(file_path[-3:])
        if file_path[-3:] == "O1s":
            if header == 'C=O':
                color = 'red'
            elif header == 'C-O':
                color = 'green'
            elif header == 'TiO2':
                color = 'blue'
            else:
                color = 'black'
        if file_path[-3:] == "C1s":
            if header == 'C-C':
                color = 'red'
            elif header == 'C=C':
                color = 'green'
            elif header == 'C-O':
                color = 'blue'
            elif header == 'C=O':
                color = 'cyan'
            elif header == 'pi-pi':
                color = 'brown'
            elif header == 'TiC':
                color = 'magenta'            
            else:
                color = 'black'        
        else:
            color = 'black'
        fit_component = df.iloc[:, i].to_numpy() + bg
        fit_area = -np.trapz(fit_component, bind)
        bg_area = -np.trapz(bg, bind)
        area_comp = fit_area-bg_area
        print(file_path,header,area_comp)
        zero = np.average(bg)
        axis.plot(bind,fit_component-zero-shift, color=color,linewidth=1.5)
    axis.plot(bind,count-zero-shift,'-',linewidth=1,color='black')
    axis.plot(bind,bg-zero-shift,color='black',linewidth=1)
    axis.plot(bind,sum_fit-zero-shift+bg,'--',color='red',linewidth=1)
    # axis.plot(bind,residual-mv_res-shift-bg,color='royalblue')
    axis.set_xlim(bind[0],bind[-1])
    # axis.set_ylim(-10000,28000)
    axis.set_xlabel('Binding energy (eV)')
    axis.set_ylabel('Intensity (arb. unit)')
    return bind, count, bg,sum_fit,residual, std_res


axis_name2 = 'ax' + '2'
axis2 = globals()[axis_name2]
bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.5-3-C1s',axis2,900,12000)
bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.1-3-C1s',axis2,900,8000)
bind, count, bg,sum_fit,residual, std_res = plot_xps('HT-2-C1s',axis2,900,4000)
bind, count, bg,sum_fit,residual, std_res = plot_xps('UT-3-C1s',axis2,900,0)



ax2.text(282, -7300, 'TiC', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
ax2.text(282, -11000, 'TiC', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
ax2.text(284.5, 26300, 'C=C', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)

ax2.annotate('C-C', xy=(285.1, 1100), xytext=(282, 2300),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('C-C', xy=(285.1, 1000-4650), xytext=(282, 1600-3300),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('C-C', xy=(285.1, 1500-4500-3600), xytext=(282, 2300-4500-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('C-C', xy=(285.1, 1500-4500-4000-4200), xytext=(282, 2300-4500-4000-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

ax2.annotate('C-O', xy=(286.35, 1100), xytext=(287, 2300),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('C-O', xy=(286.35, 400-3300), xytext=(287, 1600-3300),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('C-O', xy=(286.35, 400-3300-4250), xytext=(287, 1600-3300-4400),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('C-O', xy=(286.35, 400-3300-4300-4050), xytext=(287, 1600-3300-4400-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

# ax2.annotate('C=O', xy=(288.3, 500), xytext=(288.3, 1600),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('C=O', xy=(288.4, 300-3900), xytext=(288.4, 1600-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
# ax2.annotate('C=O', xy=(288.3, 300-4000-4000), xytext=(288.3, 1600-4200-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
# ax2.annotate('C=O', xy=(288.3, 300-4000-4000-4000), xytext=(288.3, 1600-4200-4000-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

# ax2.annotate('COOH/COOR', xy=(289, 200), xytext=(291, 2000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))
# ax2.annotate('COOH/COOR', xy=(289.3, 200-3000), xytext=(291, 2000-3100),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))
# ax2.annotate('COOH/COOR', xy=(289.3, 200-3000-3100), xytext=(291, 2000-3100-3000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))
# ax2.annotate('COOH/COOR', xy=(289.3, 200-3000-3100-3000), xytext=(291, 2000-3100-3000-3000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))

ax2.annotate('π-π*', xy=(291, 700), xytext=(291, 1500),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('π-π*', xy=(291.3, 700-3700), xytext=(291.3, 1000-3200),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('π-π*', xy=(291, 700-3700-4400), xytext=(291, 1000-3500-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax2.annotate('π-π*', xy=(291, 700-3700-4400-4000), xytext=(291, 1000-3500-4000-4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

ax2.text(292, 0+3000, 'Untreated carbon paper', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size-2)
ax2.text(292, -4000+3000, 'Heat-treated carbon paper', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size-2)
ax2.text(292, -8000+3000, 'MXene-coated (0.1mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size-2)
ax2.text(292, -12000+3000, 'MXene-coated (0.5mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size-2)

ax1.set_title('a)',y=0.995,x=0.04,pad=-14)
ax2.set_title('b) C 1s',y=0.995,x=0.1,pad=-14)
plt.tight_layout()

            # if header == 'C=O' or header == 'Ti-Cl' or header == 'Ti-C'  or header == 'TiC' or header == 'C=C' or header == 'C-Ti-F':
            #     color = 'red'
            # elif header == 'C-O' or header == 'Ti 2+' or header == 'C-O-C' or header == 'Ti-F-Ti':
            #     color = 'green'
            # elif header == 'Ti 3+' or header == 'C-C' or header == 'TiO2':
            #     color = 'blue'
            # elif header == 'C-Ti':
            #     color = 'brown'
            # else:
            #     color = 'black'