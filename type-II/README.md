This repository contains the posteriors obtained from fitting a geometric model to the dynamic radio spectrum of the M dwarf StKM 1-1262 presented in Callingham et al. (Nature, in press). The model describes the radio emission originating from a flaring loop, driven by the electron cyclotron maser instability. The posteriors can be accessed in Python as follows:
```python
import numpy as np

posteriors = np.load('posteriors.npy', allow_pickle = True).item()

theta_l = posteriors['theta_l'] # Magnetic co-latitude of flaring loop center (degrees)
phi_l = posteriors['phi_l']     # Longitude of flaring loop at t_0 (JD = 2457526.328) (degrees)
delta_l = posteriors['delta_l'] # Loop orientation relative to the meridion (degrees)
B_fp = posteriors['B_fp']       # Loop footpoint field strength (Gauss)
m = posteriors['m']             # Magnetic field gradient (Gauss/stellar radii)
L = posteriors['L']             # Loop size (stellar radii)
alpha = posteriors['alpha']     # Emission cone opening angle (degrees)
dalpha = posteriors['dalpha']   # Emission cone thickness (degrees)
flux_0 = posteriors['flux_0']   # Maximum flux density (mJy)
```
These can be easily plotted as follows:
```python
import matplotlib.pyplot as plt

plt.hist(alpha)
plt.show()
```
