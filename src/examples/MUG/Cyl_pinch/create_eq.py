#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# Create the screw pinch ("straight tokamak") equilibrium used by the MUG
# cylindrical external kink example (MUG_cyl_pinch) and print the growth rate
# predicted by Shafranov's flat-current model for the same parameters.
#
# The equilibrium has a uniform axial field, a flat current density and a flat
# density inside the plasma radius "a", both smoothed with a tanh of width "w".
# The total current is set so that the edge safety factor is "qa", making the
# field outside the plasma identical to the step (flat-current) model. The
# temperature profile provides radial force balance, dp/dr = -J_z*B_theta, with
# p = 2*n*T (Te = Ti) and a uniform temperature "Tvac" outside the plasma.
#
# Output file format (MUG_cyl_pinch "eq_file"): the first line holds the number
# of points, followed by one line per point with
#   r [m], B_theta [T], B_z [T], n [m^-3], T [eV]
#------------------------------------------------------------------------------
import argparse
import numpy as np
from scipy import integrate

mu0 = np.pi*4.E-7
proton_mass = 1.672621777E-27
elec_charge = 1.602176565E-19

parser = argparse.ArgumentParser(description="Create MUG cylindrical kink equilibrium")
parser.add_argument("--qa", type=float, default=1.6, help="Edge safety factor")
parser.add_argument("--m", type=int, default=2, help="Poloidal mode number")
parser.add_argument("--n", type=int, default=1, help="Axial mode number (k = n/R0)")
parser.add_argument("--a", type=float, default=0.1, help="Plasma radius [m]")
parser.add_argument("--b", type=float, default=0.2, help="Conducting wall radius [m]")
parser.add_argument("--R0", type=float, default=0.5, help="Equivalent major radius, L/(2*pi) [m]")
parser.add_argument("--Bz", type=float, default=0.8, help="Axial magnetic field [T]")
parser.add_argument("--n0", type=float, default=1.E20, help="Plasma density [m^-3]")
parser.add_argument("--nvac", type=float, default=0.01, help="Vacuum density fraction")
parser.add_argument("--Tvac", type=float, default=5.0, help="Vacuum temperature [eV]")
parser.add_argument("--w", type=float, default=0.004, help="Width of edge tanh [m]")
parser.add_argument("--mu_ion", type=float, default=2.0, help="Ion mass [amu]")
parser.add_argument("--npts", type=int, default=1000, help="Number of radial points")
parser.add_argument("--out_file", type=str, default="kink_eq.dat", help="Output file")
parser.add_argument("--plot", action="store_true", default=False, help="Plot profiles")
args = parser.parse_args()

# Radial grid (extends slightly past the wall)
r = np.linspace(0.0, 1.05*args.b, args.npts)
shape = 0.5*(1.0 - np.tanh((r-args.a)/args.w))

# Current density, normalized to give the desired edge safety factor
# q_a = a*B_z/(R0*B_theta(a)) with B_theta(a) = mu0*I/(2*pi*a)
Btheta_a = args.a*args.Bz/(args.R0*args.qa)
I_total = 2.0*np.pi*args.a*Btheta_a/mu0
J = shape/(2.0*np.pi*integrate.trapezoid(shape*r, r))*I_total
Ienc = np.zeros(r.shape)
Ienc[1:] = 2.0*np.pi*integrate.cumulative_trapezoid(J*r, r)
Btheta = np.zeros(r.shape)
Btheta[1:] = mu0*Ienc[1:]/(2.0*np.pi*r[1:])
Bz = args.Bz*np.ones(r.shape)

# Density
dens = args.n0*(args.nvac + (1.0-args.nvac)*shape)

# Temperature from force balance, p = p_vac + int_r^rmax J_z*B_theta dr'
pressure = np.zeros(r.shape)
pressure[:-1] = -integrate.cumulative_trapezoid((J*Btheta)[::-1], r[::-1])[::-1]
pressure += 2.0*args.n0*args.nvac*elec_charge*args.Tvac
temp = pressure/(2.0*dens*elec_charge)

# Safety factor
q = np.zeros(r.shape)
q[1:] = r[1:]*Bz[1:]/(args.R0*Btheta[1:])
q[0] = q[1]

# Growth rate from the flat-current model with a conducting wall at r=b
#   (gamma*tau_A)^2 = 2/qa^2 * (m-n*qa) * [1 - (m-n*qa)/(1-(a/b)^(2m))]
# where tau_A = R0/V_A, valid for large aspect ratio (a/R0 << 1) and low beta
rho0 = args.n0*args.mu_ion*proton_mass
V_A = args.Bz/np.sqrt(mu0*rho0)
tau_A = args.R0/V_A
alpha = np.power(args.a/args.b, 2*args.m)
delta = args.m - args.n*args.qa
gam_norm2 = 2.0*delta*(1.0 - delta/(1.0-alpha))/(args.qa**2)
print("Equilibrium parameters")
print("  a/R0      = {0:.4f}  (k*a = {1:.4f})".format(args.a/args.R0, args.n*args.a/args.R0))
print("  q(0), q(a) = {0:.4f}, {1:.4f}".format(q[0], np.interp(args.a, r, q)))
print("  I_p       = {0:.4E} [A]".format(I_total))
print("  T(0)      = {0:.4E} [eV]".format(temp[0]))
print("  beta(0)   = {0:.4E}".format(2.0*mu0*pressure[0]/args.Bz**2))
print("  V_A       = {0:.4E} [m/s]".format(V_A))
print("  tau_A     = {0:.4E} [s]".format(tau_A))
print("Flat-current model (m={0}, n={1}, b/a={2:.3f})".format(args.m, args.n, args.b/args.a))
print("  Unstable for {0:.4f} < q_a < {1:.4f}".format((args.m-1+alpha)/args.n, args.m/args.n))
if gam_norm2 > 0.0:
    gam_norm = np.sqrt(gam_norm2)
    print("  gamma*tau_A = {0:.4E}".format(gam_norm))
    print("  gamma       = {0:.4E} [1/s]".format(gam_norm/tau_A))
else:
    print("  Mode is stable")

with open(args.out_file, 'w+') as fid:
    fid.write('{0}\n'.format(args.npts))
    for i in range(args.npts):
        fid.write('{0:.12E} {1:.12E} {2:.12E} {3:.12E} {4:.12E}\n'.format(r[i], Btheta[i], Bz[i], dens[i], temp[i]))

if args.plot:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 1, sharex=True)
    ax[0].plot(r, J/1.E6)
    ax[0].set_ylabel(r'$J_z$ [MA/m$^2$]')
    ax[1].plot(r, q)
    ax[1].axhline(args.m/args.n, color='k', linestyle='--')
    ax[1].set_ylim(0.0, 4.0)
    ax[1].set_ylabel('q')
    ax[2].plot(r, dens)
    ax[2].set_ylabel(r'n [m$^{-3}$]')
    ax[3].plot(r, temp)
    ax[3].set_ylabel('T [eV]')
    ax[3].set_xlabel('r [m]')
    for ax_tmp in ax:
        ax_tmp.axvline(args.a, color='k', linestyle=':')
    plt.show()
