import h5py
import pandas as pd
import re

#convert the profiling h5 file into a CSV document for visualization in R
#Analyzing the data with Pandas would be better

dat = pd.read_hdf('starsh_gen-20-54000-1000-lfq-51S1Qw.h5', key = 'events')
dat.to_csv('starsh_gen.csv')
