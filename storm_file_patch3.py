import numpy as np
from datetime import date
from netCDF4 import Dataset
import multiprocessing
from functools import partial
from matplotlib import pyplot as plt
import time

print("Program start")

import storm_functions as storm

filename = 'D:/Robert Docs/College/NIU/GEOG 790 (SP 19)/Project/stormData/storm_track_slp'
tracked_storms = np.load(filename + '.npz', encoding='latin1', allow_pickle=True)
storms = tracked_storms['storms']

study_domain = []


np.savez('storm_track_slp_patched', storms=storms)	