"""Console script for qos."""
#%matplotlib nbagg

import matplotlib.pyplot as plt
import numpy as np
import sysid
import scipy.signal
import scipy.io


"""This file conduct output only modal analysis and finds the reference mode shapes from acceleration response data
Input: 
    y = estimated response data (row--> sensor node column --> time step
    fs = sampling rate (hz)
    sensor_pos  = sensor location (inch) for estimated data
    
    Note: When the stabilization graph is plotted, need to manually pick the poles using mouse click. Selection is
    based on user observation"""

y = np.loadtxt("estimated_acceleration.csv",delimiter=",")
fs = 512
sensor_pos=[0,60,120,180,240,300,360,420,480,540,600,660,720,780,840,900,960,1020,1080,1140,1200]

# Covariance driven Stochastic system identification
ssid = sysid.CovarianceDrivenStochasticSID(y, fs)

modes = dict()
for i, order in enumerate(range(10, 120, 2)):
    A, C, G, R0 = ssid.perform(order,25)
    modes[order] = sysid.Mode.find_modes_from_ss(A, C, fs)

sd = sysid.StabilizationDiagram()
sd.plot(modes)
f, psd = ssid.psdy(nperseg=2**12)

for i in range(2):
    freqs, Pyy = scipy.signal.csd(y[i],y[i], fs, nperseg=2 ** 12)
    sd.axes_psd.semilogy(freqs, Pyy,color=(0., 0., 0.+i, .5), lw=0.3)

plt.show()

modes = sd.picked_modes # some objects

figmodes, axes = plt.subplots(ncols=1, nrows=5, dpi=144)

ref_frequencies=[]
ref_damping=[]
ref_modeshape=[]

# Note: The reference file is to be generated using datagen to create true reference data with true natural frequnecis
# Otherwise the following comparison would not work

ref = np.load("node19analysis-ref.npz")
moderef = ref["ref_modeshapes"]
sensor_posref=ref["sensor_pos"]

for n in range(5):
    ax = axes[n]
    phi=moderef[n,:]
    #phi=np.r_[0., phi, 0.]
    mode = modes[n]
    v = np.r_[0., mode.v, 0.]
    #v = sysid.utils.modal_scale_factor(v, phi) * v
    ax.plot(sensor_pos, v.real, label='Estimated', lw=0.5, marker='x', markersize=1,c='r')
    ax.plot(sensor_posref, phi.real, label='Reference', lw=0.5, marker='x', markersize=1,c='g')
    if n == 0:
        ax.legend(bbox_to_anchor=(1, 1.20), loc=2, ncol=1)
    ref_frequencies.append(mode.f, )
    ref_modeshape.append(v.real, ) #row--> mode, col--> dof

plt.show()

## Save the reference modeshapes in .npz format
"""
np.savez('node19analysis-est',
         y=y19, fs=fs,
         ref_frequencies=ref_frequencies,
         ref_modeshapes=ref_modeshape,
         sensor_pos=sensor_pos
         )
"""

print(ref_frequencies)
