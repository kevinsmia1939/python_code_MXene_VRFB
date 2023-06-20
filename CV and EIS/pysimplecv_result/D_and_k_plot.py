import numpy as np 
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy.polynomial.polynomial as poly
import matplotlib.ticker as plticker
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Liberation Serif'
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 16 

color_lst = ["lightgreen","red","royalblue","orange","orchid"]
mark = ["o","s","^","v","D"]
fig, ax = plt.subplots(2, 2, figsize=(11,11))


def diffusion(scan,jp,alpha,conc_bulk,n):
    sqrt_scan = np.sqrt(scan)
    try: 
        jp_arr_lnfit, _ = poly.polyfit(sqrt_scan,jp,1,full=True)
        jp_arr_poly = poly.Polynomial(jp_arr_lnfit)
        jp_slope = jp_arr_lnfit[1]
        D_rev = (jp_slope/(2.69*(10**5)*n**(3/2)*conc_bulk))**2 # reversible   
        D_irr = (jp_slope/(2.99*(10**5)*n**(3/2)*(alpha**0.5)*conc_bulk))**2 # irreversible     
    except SystemError:
        pass  
    # Calculate R2
    # print(sqrt_scan)
    jp_fit = jp_arr_poly(sqrt_scan)
    residuals = jp - jp_fit
    ssr = np.sum(residuals ** 2)
    sst = np.sum((jp - np.mean(jp)) ** 2)
    r2 = (1 - (ssr / sst))
    return sqrt_scan,jp_fit,D_irr,D_rev,r2

col_idx = 1
cv_res_file = ["HT", "MX-UT-0.1",  "MX-UT-0.5"]
for i in cv_res_file:
    file = str(i) +"-neg.ods"
    df = pd.read_excel(file, engine="odf")
    jpa = df.loc[:,'Jpa'].to_numpy()
    jpc = df.loc[:,'Jpc'].to_numpy()
    scan = df.loc[:,'scanrate'].to_numpy()
    scan = (scan/1000)
    alpha = 0.5
    conc_bulk = 0.00005
    n = 1
    sqrt_scan,jp_fit,D_irr,D_rev,r2 = diffusion(scan,jpa,alpha,conc_bulk,n)
    r2 = format(r2, ".4f")
    ax[0,0].plot(sqrt_scan,jp_fit,'-',color=color_lst[col_idx])
    ax[0,0].plot(sqrt_scan,jpa,markerfacecolor=color_lst[col_idx], linestyle='', marker=mark[col_idx], markeredgecolor='black',markeredgewidth=0.8,markersize=7,label=str(i)+'   R$^2$='+str(r2))
    col_idx += 1
loc = plticker.MultipleLocator(base=1) # this locator puts ticks at regular intervals
ax[0,0].yaxis.set_major_locator(loc)    
ax[0,0].set_title('a)',y=1,x=0.05,pad=-14)
ax[0,0].legend(fontsize=12,loc='lower right',handletextpad=0)
ax[0,0].set_xlabel("v$^{1/2}$ (mV/s)$^{1/2}$")
ax[0,0].set_ylabel("J$_{pa}$ (mA/cm$^2$)")

col_idx = 1
cv_res_file = ["HT", "MX-UT-0.1",  "MX-UT-0.5"]
for i in cv_res_file:
    file = str(i) +"-neg.ods"
    df = pd.read_excel(file, engine="odf")
    jpa = df.loc[:,'Jpa'].to_numpy()
    lnjpa = np.log(jpa)
    peak_sep = df.loc[:,'ΔEₚ'].to_numpy()
    # peak_sep = -peak_sep
    jpa_fit, stats = poly.polyfit(peak_sep/2,lnjpa,1, full=True)
    jpa_poly = poly.Polynomial(jpa_fit)
    ax[0,1].plot(peak_sep*1000/2,jpa_poly(peak_sep/2),'-',color=color_lst[col_idx])
    
    
    # print((peak_sep*1000/2))
    # print(peak_sep)
    fff = np.linspace(0,max(peak_sep),10)
    # print(fff)
    # fff = np.array([124.983,142.465,156.4205])
    # ax[0,1].plot(fff*1000/2,jpa_poly(fff/2),'--',color='black',linewidth=2)
    
    #calculate R2
    residuals = lnjpa - jpa_poly(peak_sep/2)
    ssr = np.sum(residuals ** 2)
    sst = np.sum((lnjpa - np.mean(lnjpa)) ** 2)
    r_squared = (1 - (ssr / sst))
    r_squared = format(r_squared, ".4f")
    ax[0,1].plot(peak_sep*1000/2,lnjpa,markerfacecolor=color_lst[col_idx], linestyle='', marker=mark[col_idx], markeredgecolor='black',markeredgewidth=0.8,markersize=7,label=str(i)+'   R$^2$='+str(r_squared))
    col_idx += 1
    
    k_jpa = np.exp(jpa_fit[0]-np.log(0.227*96485.332*0.00005))
    # print(jpa_fit[0])
    # print('ka', k_jpa)
ax[0,1].set_title('b)',y=1,x=0.05,pad=-14)    
ax[0,1].legend(fontsize=12,loc='lower right',handletextpad=0)
ax[0,1].set_xlabel("E$_{pa}$ - E$^0$′ (mV)")
ax[0,1].set_ylabel("ln [J$_{pa}$ (mA/cm$^2$)]")
ax[0,1].set_ylim(-6.5,-5.0)
loc = plticker.MultipleLocator(base=0.5) # this locator puts ticks at regular intervals
ax[0,1].yaxis.set_major_locator(loc)
######################################################################################

col_idx = 0
cv_res_file = ["UT", "HT", "MX-UT-0.1",  "MX-UT-0.5"]
for i in cv_res_file:
    file = str(i) +"-pos.ods"
    df = pd.read_excel(file, engine="odf")
    jpa = df.loc[:,'Jpa'].to_numpy()
    scan = df.loc[:,'scanrate'].to_numpy()
    scan = scan/1000
    peak_sep = df.loc[:,'ΔEₚ'].to_numpy()
    
    jpa_fit, stats = poly.polyfit(scan**0.5,jpa,1, full=True)
    jpa_poly = poly.Polynomial(jpa_fit)
    ax[1,0].plot((scan*1000)**0.5,1000*jpa_poly(scan**0.5),'-',color=color_lst[col_idx])
    
    #calculate R2
    residuals = jpa - jpa_poly(scan**0.5)
    ssr = np.sum(residuals ** 2)
    sst = np.sum((jpa - np.mean(jpa)) ** 2)
    r_squared = (1 - (ssr / sst))
    r_squared = format(r_squared, ".4f")
    
    ax[1,0].plot((scan*1000)**0.5,1000*jpa,markerfacecolor=color_lst[col_idx], linestyle='', marker=mark[col_idx], markeredgecolor='black',markeredgewidth=0.8,markersize=7,label=str(i)+'   R$^2$='+str(r_squared))
    col_idx += 1
    D_jpa = (jpa_fit[1]/((2.99*10**5)*(0.5**0.5)*0.00005))**2 
    # print('D',D_jpa)
    
    
ax[1,0].set_title('c)',y=1,x=0.05,pad=-14)
ax[1,0].legend(fontsize=12,loc='lower right',handletextpad=0)
ax[1,0].set_xlabel("v$^{1/2}$ (mV/s)$^{1/2}$")
ax[1,0].set_ylabel("J$_{pa}$ (mA/cm$^2$)")

def reaction_rate(peak_sep,jp,alpha,conc_bulk,n):
    e_e0 = peak_sep/2
    lnjp = np.log(jp)
    try:     
        lnjp_lnfit, _ = poly.polyfit(e_e0,lnjp,1,full=True)
        lnjp_poly = poly.Polynomial(lnjp_lnfit)
        lnjpa_b = lnjp_lnfit[0] # take intercept
        k0 = np.exp(lnjpa_b-np.log(0.227*96485.332*conc_bulk))
        slope = lnjp_lnfit[1]
        # print(slope)
        
        F = 96485.332
        alpha = -slope*8.314472*298.15/F
        print(1+alpha)
    except SystemError:
        pass
    # Calculate R2
    lnjp_fit = lnjp_poly(e_e0)       
    residuals = lnjp - lnjp_fit
    ssr = np.sum(residuals ** 2)
    sst = np.sum((lnjp - np.mean(lnjp)) ** 2)
    r2 = (1 - (ssr / sst))
    # print(k0)
    return lnjp, e_e0, lnjp_fit, k0, alpha, r2

col_idx = 0
cv_res_file = ["UT", "HT", "MX-UT-0.1",  "MX-UT-0.5"]
for i in cv_res_file:
    file = str(i) +"-pos.ods"
    df = pd.read_excel(file, engine="odf")
    jp = df.loc[:,'Jpa'].to_numpy()
    peak_sep = df.loc[:,'ΔEₚ'].to_numpy()

    lnjp, e_e0, lnjp_fit, k0, alpha, r2 = reaction_rate(peak_sep,jp,alpha,conc_bulk,n)
    ax[1,1].plot(e_e0,lnjp_fit,'-',color=color_lst[col_idx])
    r2 = format(r2, ".4f")
    ax[1,1].plot(e_e0,lnjp,markerfacecolor=color_lst[col_idx], linestyle='', marker=mark[col_idx], markeredgecolor='black',markeredgewidth=0.8,markersize=7,label=str(i)+'   R$^2$='+str(r2))
    col_idx += 1
    
    
ax[1,1].set_title('d)',y=1,x=0.05,pad=-14)    
ax[1,1].legend(fontsize=12,loc='upper left',bbox_to_anchor=(0.06,1),handletextpad=0)
ax[1,1].set_xlabel("E$_{pa}$ - E$^0$′ (V)")
ax[1,1].set_ylabel("ln [J$_{pa}$ (mA/cm$^2$)]")
ax[1,1].set_ylim(-6.8,-4.7)
loc = plticker.MultipleLocator(base=0.5) # this locator puts ticks at regular intervals
ax[1,1].yaxis.set_major_locator(loc)
plt.show()
