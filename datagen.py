"""Console script for qos."""
#%matplotlib nbagg

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import sysid
from scipy import io

'''****** This file generates response data for a non-even or evenly spaced sensor network on a simply supported beam  *************

Input variables:
                pos =  Location of the sensor nodes including supports at the end
                Lt  = Total length of the beam (feet)
                Im  = Moment of inertia (in^4)
                E   =  modulus of elasticity (psi)
                m_bar   = mass density per unit length (lb.sec^2/in/in)
                T       = Data collection time (sec)
                fs      = Sampling rate (hz)
                noise   = noise level (i.e. 7*e-2 indicates 7% noise level added to the response signal
                load    = 0 for stochastic load and 1 for impact load'''

pos=np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])*12
Im= 100
E=6.58*1e7
m_bar=0.1
T=180
fs=512
noise=5*1e-2
load=0

# Change the following section only if load==1 or 2 i.e., impact load
fr=128 # impact load frequency (hz)
a=27 # unit (lb) note: 27 lb is equivalent to 120 Newton which is a reasonable impact
loc=3 # define node location of the impact
tk=2 # define exact time when impact starts (sec)

Ts=1 / fs
n=len(pos)-2

# Simulate beam
bm = sysid.utils.Beam(pos=pos, Im=Im, E=E, m_bar=m_bar)
Tmax = 1. / sysid.utils.w2f(bm.get_natural_frequency(1))
fmax = sysid.utils.w2f(bm.get_natural_frequency(n))


# Time vector
t = np.arange(0., T, Ts)

# Define loads on system
if load==0:
    ## Unmeasureable: Stochastic loads on all nodes (Process noise)
    w = np.random.normal(size=(n, t.size))*1e-1

    ## Load matrix, f
    F = w.copy()

    # Simulate response, acceleration, velocity, displacement at each node measured
    bm.set_rayleigh_damping_matrix([bm.get_natural_frequency(1), bm.get_natural_frequency(n)], [.05]*2)
    y0, _, _ = bm.simulate(t, F) # The first term indicates acceleration response under F load

    # Add measurement noise (5%)
    noise_std = y0.std()
    v = np.random.normal(size=y0.shape)* noise_std*noise
    y = y0 + v

    # Observe response signal w/ or w/o added noise
    plt.figure("Acceleration at Node 2", figsize=(10, 4))
    plt.plot(t, y[1, :], c='red')
    plt.plot(t, y0[1, :], c='orange')
    plt.xlabel("time (sec)")
    plt.ylabel("Acceleration (in/s^2")
    plt.legend(["w/ noise", "w/o noise"])
    plt.show()

elif load==1:
    ## Impact load on the system at specific nodes
    sine=a*np.sin(2*(np.pi)*fr*t)
    st=int(1/(2*fr*Ts))
    step=np.zeros((len(t)))
    step[tk*int(fs):tk*int(fs)+st]=1
    impact=sine*step

    imp=np.zeros((n,t.size))
    imp[loc]=impact # Specify which node to apply the impact load
    F=imp

    # Observe impact load
    plt.figure("Impact load", figsize=(10, 4))
    plt.plot(t, F[loc,:], c='blue')
    plt.xlabel("Load (lb)")
    plt.show()

    # Simulate response, accelerations at each floor measured
    bm.set_rayleigh_damping_matrix([bm.get_natural_frequency(1), bm.get_natural_frequency(n)], [.05]*2)
    y0, _, _ = bm.simulate(t, F)

    # Add measurement noise (5%) to the output
    noise_std = y0.std()
    v = np.random.normal(size=y0.shape)* noise_std*noise
    y = y0 + v

    # Observe response signal w/ or w/o added noise
    plt.figure("Acceleration at Node "+str(loc), figsize=(10, 4))
    plt.plot(t, y[loc, :], c='red')
    plt.plot(t, y0[loc, :], c='orange')
    plt.xlabel("time (sec)")
    plt.ylabel("Acceleration (in/s^2")
    plt.legend(["w/ noise", "w/o noise"])
    plt.show()

elif load == 2:
    ## Define Impact load with stochastic noise on the system
    sine=a*np.sin(2*(np.pi)*fr*t)
    st=int(1/(2*fr*Ts))
    step=np.zeros((len(t)))
    step[5*int(fs):5*int(fs)+st]=1
    impact=sine*step


    imp=np.zeros((n,t.size))
    imp[loc]=impact # Specify which node to apply the impact load

    ## Load matrix, f
    w = np.random.normal(size=(n, t.size))*1e-1
    F = w.copy()  # Stochastic loading
    F+= imp  # Deterministic loading to top floor

    # Observe impact load
    plt.figure("Impact load with stochastic noise", figsize=(10, 4))
    plt.plot(t, F[loc,:],c='blue')
    plt.xlabel("Load (lb)")
    plt.show()

    # Simulate response, accelerations at each floor measured
    bm.set_rayleigh_damping_matrix([bm.get_natural_frequency(1), bm.get_natural_frequency(n)], [.05] * 2)
    y0, _, _ = bm.simulate(t, F)

    # Add measurement noise (5%) to the output
    noise_std = y0.std()
    v = np.random.normal(size=y0.shape)* noise_std*noise
    y = y0 + v

    # Observe response signal w/ or w/o added noise
    plt.figure("Acceleration at Node " + str(loc), figsize=(10, 4))
    plt.plot(t, y[loc, :], c='red')
    plt.plot(t, y0[loc, :], c='orange')
    plt.xlabel("time (sec)")
    plt.ylabel("Acceleration (in/s^2")
    plt.legend(["w/ noise", "w/o noise"])
    plt.show()


plt.figure("PSD of Acceleration at last node")
for yi in [y[-1], y0[-1]]:
    freqs, Pyy = scipy.signal.csd(yi,yi, fs, nperseg=2 ** 12)
    plt.semilogy(freqs, Pyy)
    plt.xlim([0,fmax+10])

for m in range(n):
    plt.axvline(sysid.utils.w2f(bm.get_natural_frequency(m+1)), alpha=.2)

plt.ylabel('PSD')
plt.xlabel('Frequency (Hz)')
plt.show()

bm.set_rayleigh_damping_matrix([bm.get_natural_frequency(1), bm.get_natural_frequency(n)], [.05]*2) #assumed damping 5%

true_frequencies = np.array([bm.get_natural_frequency(i)/(2*np.pi) for i in range(1, n+1)])
true_damping = np.array([bm.get_rayleigh_damping_ratio(i) for i in range(1, n+1)])
true_modeshapes = np.array([bm.get_mode_shape(i) for i in range(1, n+1)])

# Save response data..............
"""
# true data in .mat format
io.savemat("full.mat",{'y':y,'fs':fs,'f':true_frequencies,'xi':true_damping,'phi':true_modeshapes,'posSensors':pos[1:-1]})

# in .csv format
np.savetxt("acc0-05.csv",
           y0,
           delimiter=",")
np.savetxt("acc-05.csv",
           y,
           delimiter=",")

"""
# in .npz format
np.savez('node19analysis-ref',
         y=y, fs=fs,
         ref_frequencies=true_frequencies,
         ref_damping=true_damping,
         ref_modeshapes=true_modeshapes,
         sensor_pos=pos
         )


print("done")