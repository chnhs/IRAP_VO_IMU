import glob
import os
import utils.utils as util
import utils.transformations as TF
import matplotlib.pyplot as plt
import numpy as np
import vo_imu_core as core
import math as m

# Find all files to run
files = glob.glob('simulated_data/*.csv')
metadata = util.PMcar()

# Run all test-cases
for idx, file in enumerate(files):
    print(f'Processing: {file}')

    # phi = m.pi/2
    # theta = m.pi/4
    # psi = m.pi/2
    # print("phi =", phi, np.rad2deg(phi))
    # print("theta  =", theta, np.rad2deg(theta))
    # print("psi =", psi, np.rad2deg(psi))
    
  
    # R = TF.Cbn_from_euler(phi, theta, psi)
    # print(np.round(R, decimals=2))
    # Euler = TF.Cbn_to_euler(R)
    # print(np.round(Euler, decimals=2))
    #exit()
    data, ground_truth = util.parse_data(file, metadata.tire_radius)
    #estimates = scheduler.run(data, metadata)
    core.VioCore(metadata.processing_mode)
    core.VioCore.run(data)
