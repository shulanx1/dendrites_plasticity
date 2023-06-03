""" Combine STA data files and save."""

import pickle
import sys
import os
import glob
wd = 'E:\\Code\\dendrites_plasticity' # working directory
sys.path.insert(1, wd)

model = 'l5'
results = wd + '\\outputs\\sta\\sta_'+model+'_data\\sta_'+model

files = glob.glob(results + '*')
cell, V, data_e, data_i, W_e, W_i = pickle.load(open(files[0], 'rb'))
for file in files[1:]:
        try:
            sta = pickle.load(open(file, 'rb'))
            V = V + sta[1]
            data_e = data_e + sta[2]
            data_i = data_i + sta[3]
            W_e = W_e + sta[4]
            W_i = W_i + sta[5]
        except:
            continue

pickle.dump([cell, V, data_e, data_i, W_e, W_i], open(wd + '\\outputs\\sta\\sta_'+model, 'wb'))
