import vamas
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.rcParams['figure.dpi'] = 100
fig, ax = plt.subplots(3, 2, figsize=(11,12))

def plot_xps(file_path,a,b,mv_res):
    df = pd.read_csv(file_path+'.csv', delimiter=',', header=1)
    size = df.shape[1]
    bind = df['corrected x'].to_numpy()
    count = df['raw_y'].to_numpy()
    sum_fit = df['sum_fit'].to_numpy()
    residual = count - sum_fit
    std_res = np.std(residual)
    bg = df['bg'].to_numpy()
    with open(file_path+'.txt', 'r') as file:
        lines = file.readlines()
        peaks = lines[13].split()[1:]
    print(file_path,peaks)
    for i in np.arange(8,size):
        header = df.columns.tolist()[i]
        if header == 'C=O' or header == 'Ti-Cl' or header == 'Ti-C'  or header == 'TiC' or header == 'C=C' or header == 'TiO2' or header == 'C-Ti-F':
            color = 'red'
        elif header == 'C-O' or header == 'Ti 2+' or header == 'C-O-C' or header == 'Ti-F-Ti':
            color = 'green'
        elif header == 'Ti 3+' or header == 'C-C':
            color = 'blue'
        elif header == 'C-Ti':
            color = 'brown'
        else:
            color = 'black'
        fit_component = df.iloc[:, i].to_numpy() + bg
        fit_area = -np.trapz(fit_component, bind)
        bg_area = -np.trapz(bg, bind)
        area_comp = fit_area-bg_area  
        zero = np.average(bg)
        
        ax[a,b].plot(bind,fit_component-zero, color=color,linewidth=0.8)
    ax[a,b].plot(bind,count-zero,'-',linewidth=0.8,color='black')
    ax[a,b].plot(bind,bg-zero,color='grey',linewidth=0.8)
    ax[a,b].plot(bind,sum_fit-zero,color='orange',linewidth=0.8)
    ax[a,b].plot(bind,residual-mv_res)
    ax[a,b].set_xlim(bind[0],bind[-1])
    ax[a,b].set_xlabel('Binding energy (eV)')
    ax[a,b].set_ylabel('Intensity (arb. unit)')
    ax[a,b].tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    return bind, count, bg,sum_fit,residual, std_res

def plot_xps_Ti(file_path,a,b,mv_res):
    df = pd.read_csv(file_path+'.csv', delimiter=',', header=1)
    size = df.shape[1]
    bind = df['corrected x'].to_numpy()
    count = df['raw_y'].to_numpy()
    sum_fit = df['sum_fit'].to_numpy()
    residual = count - sum_fit
    std_res = np.std(residual)
    bg = df['bg'].to_numpy()
    with open(file_path+'.txt', 'r') as file:
        lines = file.readlines()
        peaks = lines[11].split()[1:]
    print(peaks)
    # print(lines)
    for i in np.arange(7,size):
        header = df.columns.tolist()[i]
        if header == 'C=O' or header == 'Ti-C'  or header == 'TiC' or header == 'C=C':
            color = 'red'
        elif header == 'C-O' or header == 'Ti 2+' or header == 'C-O-C':
            color = 'green'
        elif header == 'Ti 3+' or header == 'C-C':
            color = 'blue'
        elif header == 'TiO2' or header == 'Ti 4+' or header == 'C-Ti':
            color = 'brown'
        else:
            color = 'black'
        fit_component = df.iloc[:, i].to_numpy() + bg
        fit_area = -np.trapz(fit_component, bind)
        bg_area = -np.trapz(bg, bind)
        area_comp = fit_area-bg_area  
        zero = np.average(bg)
        
        ax[a,b].plot(bind,fit_component, color=color,linewidth=0.8)
    ax[a,b].plot(bind,count,'-',linewidth=0.8,color='black')
    ax[a,b].plot(bind,bg,color='grey',linewidth=0.8)
    ax[a,b].plot(bind,sum_fit+bg,color='orange',linewidth=0.8)
    ax[a,b].plot(bind,residual-mv_res-bg)
    ax[a,b].set_xlim(bind[0],bind[-1])
    ax[a,b].set_xlabel('Binding energy (eV)')
    ax[a,b].set_ylabel('Intensity (arb. unit)')
    ax[a,b].tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
    return bind, count, bg,sum_fit,residual, std_res

excel = '../Excel sheets/' + 'MX-UT-0.1' + '.xlsx'
df = pd.read_excel(excel,sheet_name='0-wide')  
bind = df.iloc[1:,0].to_numpy()+0.4017 #Put C1s to 284.5 eV
count = df.iloc[1:,1].to_numpy()
ax[0,0].plot(bind,count,linewidth = 0.8,color='black')
ax[0,0].set_xlabel('Binding energy (eV)')
ax[0,0].set_ylabel('Intensity (arb. unit)')
ax[0,0].set_xlim(bind[0],bind[-1])
ax[0,0].set_ylim(0,450000)
ax[0,0].text(284.5+23, 400000, 'C1s', rotation=90, fontsize=12)
ax[0,0].text(450+23, 140000, 'Ti2p', rotation=90, fontsize=12)
ax[0,0].text(531+23, 120000, 'O1s', rotation=90, fontsize=12)
ax[0,0].text(200+23, 30000, 'Cl2p', rotation=90, fontsize=12)
ax[0,0].text(685+23, 120000, 'F1s', rotation=90, fontsize=12)
ax[0,0].tick_params(axis='y',which='both',left=False,right=False,labelleft=False)

# bind, count, bg,sum_fit,residual, std_res = plot_xps('UT-O1s',0,1)
bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.1-O1s',0,1,400)
ax[0,1].text(532.4, 2700, 'C-O', color='green',horizontalalignment='center',verticalalignment='center')
ax[0,1].text(529, 3100, 'TiO$_2$', color='red',horizontalalignment='center',verticalalignment='center')

bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.1-C1s',1,0,2000)
ax[1,0].text(283.6, 20000, 'C=C', color='red',horizontalalignment='center',verticalalignment='center')
ax[1,0].text(285.5, 1000, 'C-C', color='blue',horizontalalignment='center',verticalalignment='center')
ax[1,0].text(287.5, 1900, 'C-O-C', color='green',horizontalalignment='center',verticalalignment='center')
ax[1,0].text(282, 1300, 'Ti-C', color='brown',horizontalalignment='center',verticalalignment='center')
ax[1,0].text(291.1, 1100, 'π-π*', color='black',horizontalalignment='center',verticalalignment='center')

bind, count, bg,sum_fit,residual, std_res = plot_xps_Ti('MX-UT-0.1-Ti2p',1,1,-10000)
ax[1,1].text(454, 21000, 'TiC', color='red',horizontalalignment='center',verticalalignment='center')
# ax[1,1].text(456.1, 20000, 'Ti $^{2+}$', color='green',horizontalalignment='center',verticalalignment='center')
ax[1,1].annotate('Ti $^{2+}$', xy=(456.1, 19300), xytext=(458, 20000),color='green',horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))
ax[1,1].annotate('Ti $^{3+}$', xy=(457.2, 15500), xytext=(458.2, 18000),color='blue',horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))
ax[1,1].text(459, 16000, 'Ti $^{4+}$', color='brown',horizontalalignment='center',verticalalignment='center')
# ax[1,1].text(291.1, 1100, 'π-π*', color='black',horizontalalignment='center',verticalalignment='center')

bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.1-F1s',2,0,300)
ax[2,0].text(683.2, 3000, 'C-Ti-F', color='red',horizontalalignment='center',verticalalignment='center')
# ax[2,0].text(688, 850, 'Ti-F-Ti', color='green',horizontalalignment='center',verticalalignment='center')
ax[2,0].annotate('Ti-F-Ti', xy=(685.9, 530), xytext=(688, 850),color='green',horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))

bind, count, bg,sum_fit,residual, std_res = plot_xps('MX-UT-0.1-Cl2p',2,1,100)
ax[2,1].text(198.3, 550, 'Ti-Cl', color='red',horizontalalignment='center',verticalalignment='center')
ax[2,1].annotate('Ti-Cl (2p$_{1/2}$)', xy=(201.3, 340), xytext=(202, 500),horizontalalignment='center',arrowprops=dict(arrowstyle="->",facecolor='black'))

ax[0,0].set_title('a)',y=0.99,x=0.04,pad=-14)    
ax[0,1].set_title('b)',y=0.99,x=0.04,pad=-14)  
ax[1,0].set_title('c)',y=0.99,x=0.04,pad=-14)  
ax[1,1].set_title('d)',y=0.99,x=0.04,pad=-14)  
ax[2,0].set_title('e)',y=0.99,x=0.04,pad=-14)
ax[2,1].set_title('f)',y=0.99,x=0.04,pad=-14)    
plt.show()