# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:44:18 2019

@author: Joey
"""


import mat4py
import numpy as np
import matplotlib.pyplot as plt


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


for name in data.keys():
    if not name == "time":
        fig = plt.figure()
        plt.plot(data["time"], data[name])
        fig.suptitle(name)

plt.show()
        