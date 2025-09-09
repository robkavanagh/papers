This repository contains the posteriors obtained from fitting a geometric model to the dynamic radio spectrum of the M dwarf StKM 1-1262 presented in Callingham et al. (Nature, in press). The model describes the radio emission originating from a flaring loop, driven by the electron cyclotron maser instability. The posteriors can be accessed in Python as follows:
```
import numpy as np

# Load the posteriors
posteriors = np.load('posteriors.npy', allow_pickle = True).item()

theta_l = posteriors['theta_l']
phi_l = posteriors['phi_l']
delta_l = posteriors['delta_l']
B_fp = posteriors['B_fp']
m = posteriors['m']
L = posteriors['L']
alpha = posteriors['alpha']
dalpha = posteriors['dalpha']
flux_0 = posteriors['flux_0']
```
These can be easily plotted as follows:
```
import matplotlib.pyplot as plt

plt.hist(alpha)
plt.show()
```
