#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Compute the growth rate of the external kink from the history file of a
# MUG_cyl_pinch run and compare it to the flat-current model.
#
# The growth rate is fit to the kinetic energy, which grows as exp(2*gamma*t),
# over the last part of the run (after the initial transient has decayed).
#
# Usage: python plot_growth.py [--hist_file xmhd.hist] [--gamma_ref 1.0446E6]
#------------------------------------------------------------------------------
import argparse
import numpy as np
import matplotlib.pyplot as plt
from OpenFUSIONToolkit.io import histfile

parser = argparse.ArgumentParser(description="Plot growth of external kink in MUG_cyl_pinch")
parser.add_argument("--hist_file", type=str, default="xmhd.hist", help="History file")
parser.add_argument("--gamma_ref", type=float, default=1.0446E6, help="Reference growth rate [1/s] (from create_eq.py)")
parser.add_argument("--tau_A", type=float, default=4.0523E-7, help="Alfven time [s] (from create_eq.py)")
parser.add_argument("--fit_frac", type=float, default=0.5, help="Fraction of run (from end) used for fit")
parser.add_argument("--save", type=str, default=None, help="Save figure to file")
args = parser.parse_args()

hist = histfile(args.hist_file)
time = np.asarray(hist['time'])
ven = np.asarray(hist['ven'])
men = np.asarray(hist['men'])

# Fit growth rate to kinetic energy
fit_mask = time >= time[-1]*(1.0-args.fit_frac)
fit = np.polyfit(time[fit_mask], np.log(ven[fit_mask]), 1)
gamma = fit[0]/2.0
# Instantaneous growth rate
gamma_inst = 0.5*np.gradient(np.log(ven), time)

print('Growth rate (MUG)         = {0:.4E} [1/s] (gamma*tau_A = {1:.4f})'.format(gamma, gamma*args.tau_A))
print('Growth rate (flat-current) = {0:.4E} [1/s] (gamma*tau_A = {1:.4f})'.format(args.gamma_ref, args.gamma_ref*args.tau_A))
print('Relative difference        = {0:.2f}%'.format((gamma/args.gamma_ref-1.0)*100.0))

t_us = time*1.E6
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(6, 6), constrained_layout=True)
ax[0].semilogy(t_us, ven, label='Kinetic')
ax[0].semilogy(t_us, men, label='Magnetic (perturbed)')
ax[0].semilogy(t_us[fit_mask], np.exp(np.polyval(fit, time[fit_mask])), 'k--', label='Fit')
ax[0].set_ylabel('Energy [arb.]')
ax[0].legend(loc='lower right')
ax[1].plot(t_us[1:], gamma_inst[1:]*args.tau_A, label='MUG')
ax[1].axhline(args.gamma_ref*args.tau_A, color='k', linestyle='--', label='Flat-current model')
ax[1].set_ylim(0.0, 2.0*args.gamma_ref*args.tau_A)
ax[1].set_ylabel(r'$\gamma \tau_A$')
ax[1].set_xlabel(r'Time [$\mu$s]')
ax[1].legend(loc='lower right')
if args.save is not None:
    fig.savefig(args.save, dpi=150)
else:
    plt.show()
