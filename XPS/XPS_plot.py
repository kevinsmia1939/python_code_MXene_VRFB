import vamas
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.rcParams['figure.dpi'] = 100
fig, ax = plt.subplots(1, figsize=(5,10))

# shift_up = 0
# excel_file =['UT','HT','MX-UT-0.1','MX-UT-0.5','MX-HT-0.1','MX-HT-0.5']
# # excel_file =['UT','HT']
# for i in excel_file:
#     excel = 'Excel sheets/' + str(i) + '.xlsx'
#     df = pd.read_excel(excel,sheet_name='0-wide')  
#     bind = df.iloc[1:,0].to_numpy()
#     count = df.iloc[1:,1].to_numpy() + shift_up
#     plt.plot(bind,count,label=i)
#     # shift_up += 0
#     shift_up += 500000
# plt.legend(fontsize=7)
# plt.xlim(1200,0)
# plt.show()

# fig, ax = plt.subplots(1, figsize=(5,5))
# shift_up = 0
# excel_file =['UT','HT']
# for i in excel_file:
#     excel = 'Excel sheets/' + str(i) + '.xlsx'
#     df = pd.read_excel(excel,sheet_name='0-wide')  
#     bind = df.iloc[1:,0].to_numpy()
#     count = df.iloc[1:,1].to_numpy() + shift_up
#     plt.plot(bind,count,label=i)
# plt.legend(fontsize=7)
# plt.xlim(bind[0],bind[-1])

# plt.show()
# ax.cla()

plt.cla()
plt.clf()
fig, ax = plt.subplots(1, figsize=(5,5))
shift_up = 0
excel_file =['UT','HT','MX-UT-0.5']
for i in excel_file:
    excel = 'Excel sheets/' + str(i) + '.xlsx'
    df = pd.read_excel(excel,sheet_name='1-C 1s')  
    bind = df.iloc[1:,0].to_numpy()
    count = df.iloc[1:,1].to_numpy() + shift_up
    plt.plot(bind,count,label=i,linewidth=0.8)
    # shift_up += -4900
plt.legend(fontsize=7)
# plt.xlim(550,515)
plt.xlim(bind[0],bind[-1])
plt.title('Carbon 1s')
plt.show()


plt.cla()
plt.clf()
fig, ax = plt.subplots(1, figsize=(5,5))
shift_up = 0
excel_file =['UT','HT','MX-UT-0.5']
for i in excel_file:
    excel = 'Excel sheets/' + str(i) + '.xlsx'
    df = pd.read_excel(excel,sheet_name='4-O1s')  
    bind = df.iloc[1:,0].to_numpy()
    count = df.iloc[1:,1].to_numpy() + shift_up
    plt.plot(bind,count,label=i)
    # shift_up += -4900
plt.legend(fontsize=7)
# plt.xlim(550,515)
plt.xlim(bind[0],bind[-1])
plt.show()

# import vamas
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.rcParams['figure.dpi'] = 150

# data = vamas.Vamas("Casa format/UT 18 Ap location T1 in image.vms")
# # print(data)
# for block in data.blocks:
#     # print(block)
#     e_start = block.x_start
#     e_step = block.x_step
#     # print(block.corresponding_variables)
#     for var in block.corresponding_variables:
#         # print(var)
#         y_values = var.y_values
#         x_values = np.linspace( e_start,
#             (e_step * len(y_values) + e_start) - e_step,
#             num=len(y_values),
#         )
        
#         ke = x_values
#         hv = 1486.69
#         be = hv - (0) - ke
        
#         print(list(be))
#         print(y_values)
        
        
#         plt.figure()
#         # plt.xlabel(block.x_label)
#         # plt.ylabel(var.label)
#         # plt.plot(x_values, y_values)
        
#         plt.plot(be, y_values)
#         plt.xlim(be[0],be[-1])
        
        
#         plt.show()