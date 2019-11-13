import numpy as np 
import matplotlib.pyplot as plt 


def maxwell_3D(velocity, T):
    return 4*np.pi*velocity*velocity*np.exp(-velocity*velocity/(2*T)) / np.exp(1.5*np.log(2*np.pi*T))



N = 1000000
mul = 1
vx = mul*np.random.randn(N)
vy = mul*np.random.randn(N)
vz = mul*np.random.randn(N)
vel = np.sqrt(vx*vx + vy*vy + vz*vz)
Temperature = (1/3)*(vx*vx + vy*vy + vz*vz).sum() 

hist, bins = np.histogram(vel, bins=10, density=True)

bins = bins[:-1] + np.diff(bins) / 2