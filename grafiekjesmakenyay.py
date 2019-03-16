# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:11:10 2019

@author: Thomas
"""





def plot_real_data(tstart, tend, value, data):
    def knippen(tmin_s,tsec_s,tmin_e,tsec_e,lst,tlst):  
        starttime = 60*tmin_s + tsec_s
        endtime= 60*tmin_e + tsec_e
        newlist=lst[starttime*10:endtime*10+1]
        newtlist=tlst[starttime*10:endtime*10+1]
        return newlist, newtlist    
    
    #fig, ax = plt.subplots(1,1)
    
    phugoidstart = tstart
    phugoidend = tend
    
    delst,tlst = knippen(phugoidstart[0],phugoidstart[1],phugoidend[0],phugoidend[1],data[value],data["time"])
    
    #ax.plot(tlst,delst, label="test data")
    #fig.suptitle(value)
    
    return tlst, delst