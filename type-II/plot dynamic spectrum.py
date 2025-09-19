import numpy as np
import matplotlib.pyplot as plt
import json
import pickle
import astropy.units as u

rad_per_deg = np.pi / 180

##########################################
# Setup
##########################################

# Rotation period (days)
P = 1.24

# Line of sight
x, y, z = np.eye(3)
x = np.expand_dims(x, axis = (0, 1))
y = np.expand_dims(y, axis = (0, 1))
z = np.expand_dims(z, axis = (0, 1))


##########################################
# Dynamic spectrum
##########################################

# Dynamic spectrum image around the burst
data = pickle.load(open('burst.pkl', 'rb'))

time = data['Time (days)']
time = np.expand_dims(time, axis = (1, 2))

freq = data['Frequency (MHz)']
freq = np.expand_dims(freq, axis = (0, 2))

flux_obs = data['Stokes V (mJy)'] - 50
flux_noise = 39

# Rotation phase
phi_rot = 2 * np.pi * time / P

# Field strengths
B = freq / 2.8
B_min = freq.min() / 2.8


##########################################
# Model
##########################################

def model(params):

	theta_l, phi_l, delta_l, B_fp, m, L, alpha, dalpha, flux_0 = params

	B_fp_max_max = B_min + m * L

	# Encourage solution 
	if B_fp_max_max < B_fp: return -1e50 * abs(B_fp_max_max - B_fp)

	# Longitude of loop as function of time
	phi = phi_l + phi_rot

	# Repeated trig terms
	sin_theta_l = np.sin(theta_l)
	cos_theta_l = np.cos(theta_l)
	sin_delta_l = np.sin(delta_l)
	cos_delta_l = np.cos(delta_l)
	sin_phi = np.sin(phi)
	cos_phi = np.cos(phi)

	# Vectors
	n_l = sin_phi * y + cos_phi * x
	y_l = cos_phi * y - sin_phi * x
	x_l = cos_theta_l * z + sin_theta_l * n_l
	z_l = sin_theta_l * z - cos_theta_l * n_l
	d_l = sin_delta_l * y_l + cos_delta_l * z_l

	# Compute flux from loop vector orientation
	sin_theta = (1 / (2 * L)) * ((B_fp - B) / m + 1) ** 2 - L / 2 - 1 / (2 * L)
	cos_theta = (1 - sin_theta ** 2) ** 0.5
	c_l = cos_theta * x_l - sin_theta * d_l
	gamma = np.arccos(c_l[:, :, 0])
	return flux_0 * np.exp(- 0.5 * (((gamma - alpha) / dalpha)) ** 2)


##########################################
# Plot
##########################################

##########################################
# Figure setup
##########################################

fig_width = 6
fig_height = 4.5

left = 0.1
bottom = 0.1
right = 0.1
top = 0.03 * fig_width / fig_height

# hspace = 0.01
# vspace = hspace * fig_width / fig_height

hspace = 0.03
ax_cbar_width = 0.05

ax_width = 1 - left - right - hspace - ax_cbar_width
ax_height = 1 - bottom - top

fig = plt.figure(figsize = (fig_width, fig_height))
ax = fig.add_axes([left, bottom, ax_width, ax_height])
ax_cbar = fig.add_axes([1 - right - ax_cbar_width, bottom, ax_cbar_width, ax_height])

ax.set_xlabel('Time (minutes)')
ax.set_ylabel('Frequency (MHz)')
ax.set_xlim(-1.5, 1.5)

##########################################
# Plot
##########################################

data = json.load(open('results.json'))
logz = data['logz']
logl = data['maximum_likelihood']['logl']
params = data['maximum_likelihood']['point']

flux_model = model(params)

time = (time.flatten() - 0.3772) * u.day.to(u.min)

c = ax.pcolormesh(time.flatten(), freq.flatten(), flux_model.T, vmin = 0, vmax = 100, cmap = 'Blues_r', rasterized = True)
plt.colorbar(c, cax = ax_cbar, label = 'Flux density (mJy)')

plt.show()