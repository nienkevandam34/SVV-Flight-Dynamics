# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:11:10 2019

@author: Thomas
"""

import mat4py
import numpy as np
import matplotlib.pyplot as plt

def knippen(tmin_s,tsec_s,tmin_e,tsec_e,lst,tlst):  
    starttime = 60*tmin_s + tsec_s
    endtime= 60*tmin_e + tsec_e
    newlist=lst[starttime*10:endtime*10+1]
    newtlist=tlst[starttime*10:endtime*10+1]
    return newlist, newtlist

raw_data = mat4py.loadmat('reference_data.mat')["flightdata"]

data = {}
unit = {}

for name in raw_data.keys():
    print(name)
    
    new_data = []
    
    if not name == "time":
        # name.shape = (N,1)
        for old_data in raw_data[name]["data"]:
            new_data.append(old_data[0])
    else:
        # time.shape = (1,N)
        data[name] = np.array(raw_data[name]["data"])
        continue
    
    data[name] = np.array(new_data)
    unit[name] = raw_data[name]["units"]


# =============================================================================
# for name in data.keys():
#     if not name == "time":
#         fig = plt.figure()
#         plt.plot(data["time"], data[name])
#         fig.suptitle(name)
# =============================================================================


phugoidstart = 53,0
phugoidend = 58,0

delst,tlst = knippen(phugoidstart[0],phugoidstart[1],phugoidend[0],phugoidend[1],data["Dadc1_mach"],data["time"])


# =============================================================================
# for name in data.keys():
#     if not name = "time"
# =============================================================================

plt.plot(tlst,delst)

plt.show()