#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Plot the structure of the external kink from a MUG_cyl_pinch plot run
# (plot_run=T), using the plot files (oft_xdmf.*.h5) and the magnetic probe
# signals at the plasma edge (b_probes.hist).
#
# Left: perturbed radial magnetic field (color) and perturbed velocity inside
# the plasma (arrows) in the z=0 plane. Right: perturbed radial magnetic field
# around the plasma edge, which shows the m=2 structure of the mode.
#
# Usage: python plot_mode.py [--a 0.1] [--save MUG_cyl_kink_ex-Mode.png]
#------------------------------------------------------------------------------
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from OpenFUSIONToolkit.io import XDMF_plot_file, histfile

parser = argparse.ArgumentParser(description="Plot mode structure of external kink in MUG_cyl_pinch")
parser.add_argument("--a", type=float, default=0.1, help="Plasma radius [m]")
parser.add_argument("--rv", type=float, default=0.125, help="Outer radius of fluid vacuum region [m]")
parser.add_argument("--save", type=str, default=None, help="Save figure to file")
args = parser.parse_args()

# Collect perturbed fields at the last plot time in the z=0 plane from all processes
pts = []
Bvals = []
Vvals = []
for filename in sorted(glob.glob('oft_xdmf.*.h5')):
    mesh = XDMF_plot_file(filename)['mug']['vmesh']
    mask = np.abs(mesh.r[:, 2]) < 1.E-8
    pts.append(mesh.r[mask, :2])
    Bvals.append(mesh.get_field('B', timestep=len(mesh.times)-1)[mask, :])
    Vvals.append(mesh.get_field('V', timestep=len(mesh.times)-1)[mask, :])
    plot_time = mesh.times[-1]
pts, inds = np.unique(np.round(np.vstack(pts), 10), axis=0, return_index=True)
Bvals = np.vstack(Bvals)[inds, :]
Vvals = np.vstack(Vvals)[inds, :]
theta = np.arctan2(pts[:, 1], pts[:, 0])
Br = Bvals[:, 0]*np.cos(theta) + Bvals[:, 1]*np.sin(theta)
Br /= np.abs(Br).max()

# Perturbed radial field from probes around the plasma edge
probes = histfile('b_probes.hist')
Bprobe = np.asarray(probes['B'])[-1, :].reshape((-1, 3))
theta_p = np.linspace(0.0, 2.0*np.pi, Bprobe.shape[0], endpoint=False)
Br_p = Bprobe[:, 0]*np.cos(theta_p) + Bprobe[:, 1]*np.sin(theta_p)
Br_p /= np.abs(Br_p).max()
modes = np.abs(np.fft.rfft(Br_p))
print('Fraction of edge B_r power in m=2: {0:.4f}'.format(modes[2]**2/np.sum(modes[1:]**2)))

fig, ax = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
tri = mtri.Triangulation(pts[:, 0], pts[:, 1])
clf = ax[0].tricontourf(tri, Br, levels=np.linspace(-1.0, 1.0, 21), cmap='RdBu_r')
fig.colorbar(clf, ax=ax[0], label=r'$\delta B_r$ [norm.]')
inner = np.sqrt(pts[:, 0]**2 + pts[:, 1]**2) < args.a
vscale = np.sqrt(np.sum(Vvals[inner, :2]**2, axis=1)).max()
skip = max(1, int(inner.sum()/300))
ax[0].quiver(pts[inner, 0][::skip], pts[inner, 1][::skip], Vvals[inner, 0][::skip]/vscale, Vvals[inner, 1][::skip]/vscale,
             scale=20, width=0.004)
for rad, style in ((args.a, 'k--'), (args.rv, 'k:')):
    ax[0].plot(rad*np.cos(theta_p), rad*np.sin(theta_p), style, linewidth=1)
ax[0].set_aspect('equal')
ax[0].set_xlabel('x [m]')
ax[0].set_ylabel('y [m]')
ax[0].set_title(r'$z = 0$, t = {0:.1f} $\mu$s'.format(plot_time*1.E6))
ax[1].plot(np.degrees(theta_p), Br_p, label='MUG')
phase = np.arctan2(-np.sum(Br_p*np.sin(2*theta_p)), np.sum(Br_p*np.cos(2*theta_p)))
ax[1].plot(np.degrees(theta_p), np.cos(2.0*theta_p + phase), 'k--', label=r'$\cos(2\theta + \phi)$')
ax[1].set_xlim(0.0, 360.0)
ax[1].set_xticks(np.arange(0, 361, 90))
ax[1].set_xlabel(r'$\theta$ [deg]')
ax[1].set_ylabel(r'$\delta B_r(r=a)$ [norm.]')
ax[1].set_ylim(-1.1, 1.4)
ax[1].legend(loc='upper center', ncol=2)
ax[1].grid(True)
if args.save is not None:
    fig.savefig(args.save, dpi=150)
else:
    plt.show()
