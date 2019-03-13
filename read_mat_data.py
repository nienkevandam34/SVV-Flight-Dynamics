# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:44:18 2019

@author: Joey

Use this file to read the data in a .mat flight data file.
"""

def read_mat(filename):
    """
    Read a .mat file formatted in the same way as the reference data. 
    
    INPUT:
        filename (str): name of the .mat file to be read (including .mat)
        
    OUTPUT:
        data, unit, description, keys
        
    Data, unit and description are dicts, keys is a list. All dictionaries 
    can be accessed with the same keys. Data contains the data as a numpy 
    array.
    """
    
    import mat4py
    import numpy as np
    
    print("Reading {} ...".format(filename))
    raw_data = mat4py.loadmat(filename)["flightdata"]
    
    print("Sorting data ...")
    data = {}
    unit = {}
    description = {}
    keys = []
    
    for name in raw_data.keys():
        
        new_data = []
        
        if not name == "time":
            # name.shape = (N,1)
            for old_data in raw_data[name]["data"]:
                new_data.append(old_data[0])
                
            data[name] = np.array(new_data)
            
        else:
            # time.shape = (1,N)
            data[name] = np.array(raw_data[name]["data"])
        
        unit[name] = raw_data[name]["units"]
        description[name] = raw_data[name]["description"]
        keys.append(name)
    
    return data, unit, description, keys

























        