import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 100

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16 
fig, ax = plt.subplots(figsize=(15,10))

label_size = 13

ax1 = plt.subplot(131)
ax1.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax1.set_xlabel('Binding energy (eV)')
ax1.set_ylabel('Intensity (arb. unit)')

# ax2 = plt.subplot(222)
# ax2.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
# ax2.set_xlabel('Binding energy (eV)')
# ax2.set_ylabel('Intensity (arb. unit)')

ax3 = plt.subplot(132)
ax3.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax3.set_xlabel('Binding energy (eV)')
ax3.set_ylabel('Intensity (arb. unit)')

ax4 = plt.subplot(133)
ax4.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax4.set_xlabel('Binding energy (eV)')
ax4.set_ylabel('Intensity (arb. unit)')

xps_lst = ['UT','HT','MX-UT-0.1','MX-UT-0.5']
color_lst = ['green','red','royalblue','orange']
shift_idx = 0
shift=[-44000,-40000,-20000,0]
shift_l = 0

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
        peaks = lines[13].split()[1:]
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
        print(header)
        # print(file_path[-3:])
        if file_path[-3:] == "O1s":
            if header == 'C-O':
                color = 'red'
            elif header == 'C-Ti-(OH)x':
                color = 'green'
            elif header == 'TiO2':
                color = 'blue'
            elif header == 'C=O':
                color = 'magenta'
            else:
                color = 'black'
        elif file_path[-3:] == "C1s":
            if header == 'C-C':
                color = 'red'
            elif header == 'C=C':
                color = 'green'
            elif header == 'C-O':
                color = 'blue'
            elif header == 'C=O':
                color = 'yellow'
            elif header == 'O-C=O':
                color = 'brown'
            else:
                color = 'black'
        elif file_path[-4:] == "Cl2p":
            if header == 'Ti-Cl':
                color = 'green'
            else:
                color = 'blue'  
        else:
            color = 'black'
        fit_component = df.iloc[:, i].to_numpy() + bg
        fit_area = -np.trapz(fit_component, bind)
        bg_area = -np.trapz(bg, bind)
        area_comp = fit_area-bg_area
        print(file_path,header,area_comp)
        zero = np.average(bg)
        # print(color)
        axis.plot(bind,fit_component-zero-shift, color=color,linewidth=1)
    axis.plot(bind,count-zero-shift,'-',linewidth=0.8,color='black')
    axis.plot(bind,bg-zero-shift,color='black',linewidth=1)
    axis.plot(bind,sum_fit-zero-shift,'--',color='red',linewidth=1)
    # axis.plot(bind,residual-mv_res-shift,color='royalblue')
    axis.set_xlim(bind[0],bind[-1])
    # axis.set_ylim(-7000,500)
    axis.set_xlabel('Binding energy (eV)')
    axis.set_ylabel('Intensity (arb. unit)')
    return bind, count, bg,sum_fit,residual, std_res

def plot_xps_Ti(file_path,axis,mv_res,shift):
    df = pd.read_csv(file_path+'.csv', delimiter=',', header=1)
    size = df.shape[1]
    bind = df['corrected x'].to_numpy()
    count = df['raw_y'].to_numpy()
    sum_fit = df['sum_fit'].to_numpy()
    residual = count - sum_fit
    std_res = np.std(residual)
    bg = df['bg'].to_numpy()
    print("##############################################")
    with open(file_path+'.txt', 'r') as file:
        lines = file.readlines()
        peaks = lines[11].split()[1:]
    print(peaks)
    print(file_path)
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
        print(header)
        
        if file_path[-3:] == "F1s":
            if header == 'Ti-F':
                color = 'blue'
            elif header == 'TiF2':
                color = 'green'
            else:
                color = 'black'
        elif file_path[-3:] == "Ti2p":
            if header == 'C-C':
                color = 'red'
            elif header == 'C=C':
                color = 'green'
            elif header == 'C-O':
                color = 'blue'
            elif header == 'C=O':
                color = 'yellow'
            elif header == 'O-C=O':
                color = 'brown'
            else:
                color = 'black'        
        else:
            color = 'black'
        fit_component = df.iloc[:, i].to_numpy() + bg
        fit_area = -np.trapz(fit_component, bind)
        bg_area = -np.trapz(bg, bind)
        area_comp = fit_area-bg_area  
        zero = np.average(bg)
        
        axis.plot(bind,fit_component-shift, color=color,linewidth=1)
    axis.plot(bind,count-shift,'-',linewidth=1,color='black')
    axis.plot(bind,bg-shift,color='black',linewidth=1)
    axis.plot(bind,sum_fit+bg-shift,'--',color='red',linewidth=1)
    # axis.plot(bind,residual-mv_res-bg-shift,color='royalblue')
    axis.set_xlim(bind[0],bind[-1])
    # axis.set_ylim(bind[0],bind[-1])
    axis.set_xlabel('Binding energy (eV)')
    axis.set_ylabel('Intensity (arb. unit)')
    return bind, count, bg,sum_fit,residual, std_res

def plot_xps_Ti2(file_path,axis,mv_res,shift):
    df = pd.read_csv(file_path+'.csv', delimiter=',', header=1)
    size = df.shape[1]
    bind = df['corrected x'].to_numpy()
    count = df['raw_y'].to_numpy()
    sum_fit = df['sum_fit'].to_numpy()
    residual = count - sum_fit
    std_res = np.std(residual)
    bg = df['bg'].to_numpy()
    print("##############################################")
    with open(file_path+'.txt', 'r') as file:
        lines = file.readlines()
        peaks = lines[11].split()[1:]
    print(peaks)
    print(file_path)
    # print(df.columns.tolist())
    df = df.drop('raw_x', axis=1)
    df = df.drop('raw_y', axis=1)
    df = df.drop('corrected x', axis=1)
    df = df.drop('data-bg', axis=1)
    df = df.drop('sum_components', axis=1)
    df = df.drop('bg', axis=1)
    df = df.drop('sum_fit', axis=1)
    # print(df.columns.tolist())
    try:
        df = df.drop('bg_shirley_', axis=1)
    except:
        pass
    size = df.shape[1]
    for i in np.arange(0,size):
        header = str(df.columns.tolist()[i])
        print(header)
        if file_path[-4:] == "Ti2p":
            if header == 'TiC':
                color = 'red'
            elif header == 'Ti 2+':
                color = 'green'
            elif header == 'Ti 3+':
                color = 'blue'
            elif header == 'Ti 4+':
                color = 'magenta'
            else:
                color = 'black'        
        else:
            color = 'black'
        fit_component = df.iloc[:, i].to_numpy() + bg
        fit_area = -np.trapz(fit_component, bind)
        bg_area = -np.trapz(bg, bind)
        area_comp = fit_area-bg_area  
        zero = np.average(bg)
        
        axis.plot(bind,fit_component-shift, color=color,linewidth=1)
    axis.plot(bind,count-shift,'-',linewidth=1,color='black')
    axis.plot(bind,bg-shift,color='black',linewidth=1)
    axis.plot(bind,sum_fit+bg-shift,'--',color='red',linewidth=1)
    # axis.plot(bind,residual-mv_res-bg-shift,color='royalblue')
    axis.set_xlim(bind[0],bind[-1])
    # axis.set_ylim(bind[0],bind[-1])
    axis.set_xlabel('Binding energy (eV)')
    axis.set_ylabel('Intensity (arb. unit)')
    return bind, count, bg,sum_fit,residual, std_res

# ax1 = plt.subplot(422)
axis_name1 = 'ax' + '1'
axis1 = globals()[axis_name1]
os1_sh = 6000
bind, count, bg,sum_fit,residual, std_res = plot_xps('UT-O1s',axis1,600,0)
bind, count, bg,sum_fit,residual, std_res = plot_xps('HT-O1s',axis1,600,os1_sh)
bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.1-1-O1s',axis1,600,2*os1_sh)
bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.5-1-O1s',axis1,600,3*os1_sh)
axis1.text(532.8, 5400, 'C-O', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)

axis1.text(533, 5500-os1_sh-1200, 'C-O', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
axis1.annotate('C=O', xy=(531.3, 5500-os1_sh-4200), xytext=(531.3, 5500-os1_sh-2000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

axis1.annotate('C-O', xy=(532.8, 5500-2*os1_sh-4100), xytext=(532.8, 5500-2*os1_sh-1300),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
axis1.annotate('C-Ti-O$_x$ /\n C-Ti-(OH)$_x$', xy=(531.5, 5500-2*os1_sh-4000), xytext=(531.5, 5500-2*os1_sh-2200),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size-3)
axis1.annotate('Ti-O', xy=(529.2, 5500-2*os1_sh-3000), xytext=(527, 5500-2*os1_sh-2200),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

axis1.annotate('C-O', xy=(532.8, 5500-3*os1_sh-4200), xytext=(532.8, 5500-3*os1_sh-1300),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
axis1.annotate('C-Ti-O$_x$ /\n C-Ti-(OH)$_x$', xy=(531.7, 5500-3*os1_sh-4000), xytext=(531.7, 5500-3*os1_sh-2200),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size-3)
axis1.annotate('Ti-O', xy=(529.2, 5500-3*os1_sh-3000), xytext=(527, 5500-3*os1_sh-2200),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

axis1.text(536, 2000, 'Untreated \n carbon paper', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
axis1.text(536, 2000-os1_sh, 'Heat-treated \n carbon paper', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
axis1.text(535.5, 2000-2*os1_sh+1000, 'MXene-coated carbon \n paper (0.1mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
axis1.text(535.5, 2000-3*os1_sh+1000, 'MXene-coated carbon \n paper (0.5mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)


# axis_name2 = 'ax' + '2'
# axis2 = globals()[axis_name2]
# bind, count, bg,sum_fit,residual, std_res = plot_xps_Ti2('MX-UT-0.1-3-Ti2p',axis2,-9000,0)
# bind, count, bg,sum_fit,residual, std_res = plot_xps_Ti2('MX-UT-0.5-4-Ti2p',axis2,-9000,20000)
# axis2.text(468, 18000, 'MX-0.1', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
# axis2.text(468, -100, 'MX-0.5', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)

# axis2.annotate('TiC', xy=(455, -300+18000), xytext=(452, 1000+19000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
# axis2.annotate('Ti$^{2+}$', xy=(455.5, -300+19000), xytext=(453, 4000+17000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
# axis2.annotate('Ti$^{3+}$', xy=(457, -1500+16800), xytext=(458, 4000+15000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)

# axis2.annotate('TiC', xy=(455, -300), xytext=(452, 1000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
# axis2.annotate('Ti$^{2+}$', xy=(455.5, -300), xytext=(453, 4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
# axis2.annotate('Ti$^{3+}$', xy=(457, -1500), xytext=(458, 4000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
# axis2.text(461.5, 19000, 'Ti 2p$_{1/2}$', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
# axis2.text(461.5, 2000, 'Ti 2p$_{1/2}$', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
# axis2.text(455.5, 23000, 'Ti 2p$_{3/2}$', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
# axis2.text(455.5, 8000, 'Ti 2p$_{3/2}$', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
# axis2.set_ylim(-10000,24800)

axis_name3 = 'ax' + '3'
axis3 = globals()[axis_name3]
bind, count, bg,sum_fit,residual, std_res = plot_xps_Ti('MX-UT-0.1-F1s',axis3,-6800,0)
ax3.annotate('Ti-F$_2$', xy=(684.2, 7700), xytext=(683, 9500),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax3.annotate('Ti-F', xy=(686, 9000), xytext=(687, 9500),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
bind, count, bg,sum_fit,residual, std_res = plot_xps_Ti('MX-UT-0.5-F1s',axis3,-7500,7000)
ax3.annotate('Ti-F$_2$', xy=(684, 1500), xytext=(683, 3000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax3.annotate('Ti-F', xy=(686, 3000), xytext=(687, 3500),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax3.text(689,9000, 'MXene-coated carbon \n paper (0.1mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
ax3.text(689, 3000, 'MXene-coated carbon \n paper (0.5mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)



axis_name4 = 'ax' + '4'
axis4 = globals()[axis_name4]
bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.1-Cl2p',axis4,150,0)
ax4.annotate('Ti-Cl (2p$_{3/2}$)', xy=(198.8, 300), xytext=(197,400),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax4.annotate('Ti-Cl (2p$_{1/2}$)', xy=(201.1, 300), xytext=(203.5,300),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.5-Cl2p',axis4,150,1000)
ax4.annotate('Ti-Cl (2p$_{3/2}$)', xy=(198.8, 300-1000), xytext=(197,400-1000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax4.annotate('Ti-Cl (2p$_{1/2}$)', xy=(201.4, 300-1000), xytext=(203.5,300-1000),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'), fontsize=label_size)
ax4.text(205,500, 'MXene-coated carbon \n paper (0.1mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)
ax4.text(205, -500, 'MXene-coated carbon \n paper (0.5mg/cm$^2$)', color='black',horizontalalignment='center',verticalalignment='center', fontsize=label_size)


plt.tight_layout()

ax1.set_title('O 1s',y=0.995,x=0.1,pad=-14)
# ax2.set_title('b) Ti 2p',y=0.995,x=0.1,pad=-14)
ax3.set_title('F 1s',y=0.995,x=0.1,pad=-14)
ax4.set_title('Cl 2p',y=0.995,x=0.1,pad=-14)