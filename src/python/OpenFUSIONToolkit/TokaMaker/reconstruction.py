'''! Python interface for TokaMaker Grad-Shafranov functionality

@authors Chris Hansen
@date May 2023
@ingroup doxy_oft_python
'''
import numpy as np
from ..util import *

class tokamaker_recon_settings_struct(c_struct):
    r'''! TokaMaker reconstruction settings structure

     - `pm` Print 'performance' information (eg. iteration count) during run?
     - `free_boundary` Perform free-boundary calculation?
     - `has_plasma` Include plasma effects in calculation, vacuum otherwise?
     - `limited_only` Do not search for X-points when determining LCFS?
     - `maxits` Maximum NL iteration count for G-S solver
     - `mode` Parallel current source formulation used (0 -> define \f$F'\f$, 1 -> define \f$F*F'\f$)
     - `urf` Under-relaxation factor for NL fixed-point iteration
     - `nl_tol` Convergence tolerance for NL solver
     - `rmin` Minimum magnetic axis major radius, used to catch 'lost' equilibria
     - `lim_zmax` Maximum vertical range for limiter points, can be used to exclude complex diverter regions
     - `limiter_file` File containing additional limiter points not included in mesh (default: 'none')
    '''
    _fields_ = [("fitI", c_bool),
                ("fitP", c_bool),
                ("fitPnorm", c_bool),
                ("fitAlam", c_bool),
                ("fitR0", c_bool),
                ("fitV0", c_bool),
                ("fitF0", c_bool),
                ("fixedCentering", c_bool),
                ("pm", c_bool)]


def tokamaker_recon_default_settings():
    '''! Initialize reconstruction settings object with default values

    @result tokamaker_recon_settings_struct object
    '''
    settings = tokamaker_recon_settings_struct()
    settings.fitI = False
    settings.fitP = False
    settings.fitPnorm = True
    settings.fitAlam = True
    settings.fitR0 = False
    settings.fitV0 = False
    settings.fitF0 = False
    settings.fixedCentering = False
    settings.pm = False
    return settings

## @cond
tokamaker_recon_run = ctypes_subroutine(oftpy_lib.tokamaker_recon_run,
    [c_bool, ctypes.POINTER(tokamaker_recon_settings_struct), c_void_p, c_int_ptr])
## @endcond

Mirnov_con_id = 1
Ip_con_id = 2
fluxLoop_con_id = 7
dFlux_con_id = 8
Press_con_id = 9
q_con_id = 10
saddle_con_id = 11


class Mirnov_con:
    def __init__(self, pt=None, phi=0., norm=None, val=None, err=None):
        self.pt = pt
        self.phi = phi
        self.norm = norm
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        self.phi = float(values[2])
        values = file.readline().split()
        self.norm = (float(values[0]), float(values[1]), float(values[2]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(Mirnov_con_id))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.pt[0], self.pt[1], self.phi))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.norm[0], self.norm[1], self.norm[2]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class Ip_con:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(Ip_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class fluxLoop_con:
    def __init__(self, pt=None, val=None, err=None):
        self.pt = pt
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(fluxLoop_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt[0], self.pt[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class dFlux_con:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(dFlux_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class Press_con:
    def __init__(self, pt=None, val=None, err=None):
        self.pt = pt
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(Press_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt[0], self.pt[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class q_con:
    def __init__(self, type=None, val=None, err=None, loc=0.):
        self.type = type
        self.val = val
        self.err = err
        self.loc = loc

    def read(self, file):
        values = file.readline().split()
        self.type = int(values[0])
        self.loc = float(values[1])
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(q_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.type, self.loc))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class saddle_con:
    def __init__(self, pt1=None, pt2=None, width=None, val=None, err=None):
        self.pt1 = pt1
        self.pt2 = pt2
        self.width = width
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt1 = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.pt2 = (float(values[0]), float(values[1]))
        value = file.readline()
        self.width = float(value)
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(saddle_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt1[0], self.pt1[1]))
        file.write(' {0:E} {1:E}\n'.format(self.pt2[0], self.pt2[1]))
        file.write(' {0:E}\n'.format(self.width))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


con_map = {
    Mirnov_con_id: Mirnov_con,
    Ip_con_id: Ip_con,
    fluxLoop_con_id: fluxLoop_con,
    dFlux_con_id: dFlux_con,
    Press_con_id: Press_con,
    q_con_id: q_con,
    saddle_con_id: saddle_con
}


class reconstruction():
    def __init__(self,gs_obj,filename=None):
        ## Grad-Shafranov object for reconstruction
        self._gs_obj = gs_obj
        ## Reconstruction specific settings object
        self.settings = tokamaker_recon_default_settings()
        ## Plasma current constraint
        self._Ip_con = None
        ## Diamagnetic flux constraint
        self._Dflux_con = None
        ## Flux loop constraints
        self._flux_loops = []
        ## Mirnov constraints
        self._mirnovs = []
        ## Saddle loop constraints
        self._saddles = []
        ## Pressure constraints
        self._pressure_cons = []
        ## Coil current targets
        self._coil_targets = None
        ## Coil current error weights
        self._coil_wts = None
        #
        if filename is not None:
            self.read_fit_in(filename)
    
    def __del__(self):
        self._gs_obj = None
        self.settings = None
        self._Ip_con = None
        self._Dflux_con = None
        self._flux_loops = []
        self._mirnovs = []
        self._saddles = []
        self._pressure_cons = []
    
    def set_Ip(self,Ip,err):
        '''! Add a constraint on the toroidal plasma current

        @param val Current constraint [A]
        @param err Uncertainty in constraint [A]
        '''
        if Ip < 0.0:
            raise ValueError('Toroidal current must be >= 0')
        self._Ip_con = Ip_con(val=Ip, err=err)

    def set_DFlux(self,DFlux,err):
        '''! Add a constraint on the diamagnetic flux

        @param val Flux constraint [Wb]
        @param err Uncertainty in constraint [Wb]
        '''
        self._Dflux_con = dFlux_con(val=DFlux, err=err)
    
    def set_coil_currents(self, coil_currents, err=None):
        if coil_currents is None:
            self._coil_targets = None
        else:
            if coil_currents.shape[0] != self._gs_obj.ncoils:
                raise ValueError('Incorrect size for "coil_currents", must be [ncoils,]')
        if err is None:
            self._coil_wts = None
        else:
            if err.shape[0] != self._gs_obj.ncoils:
                raise ValueError('Incorrect size for "err", must be [ncoils,]')
            self._coil_wts = 1.0/abs(err)

    def add_flux_loop(self,loc,val,err):
        '''! Add a constraint for the flux measured by a full poloidal flux loop

        @param loc (R,Z) position of constraint [m]
        @param val Flux constraint [Wb]
        @param err Uncertainty in constraint [Wb]
        '''
        self._flux_loops.append(fluxLoop_con(pt=loc, val=val, err=err))

    def add_Mirnov(self,loc,norm,val,err):
        r'''! Add a constraint for the magnetic field at a point

        @param loc (R,Z) position of constraint [m]
        @param norm (R,\f$ \phi \f$,Z) direction of magnetic field
        @param val Field constraint [T]
        @param err Uncertainty in constraint [T]
        '''
        self._mirnovs.append(Mirnov_con(pt=loc, norm=norm, val=val, err=err))
    
    def add_saddle(self,p1,p2,width,val,err):
        '''! Add a constraint for the flux measured by a saddle loop

        @param p1 (R,Z) position of first toroidal leg of saddle [m]
        @param p2 (R,Z) position of second toroidal leg of saddle [m]
        @param val Flux constraint [Wb]
        @param err Uncertainty in constraint [Wb]
        '''
        self._saddles.append(saddle_con(p1=p1, p2=p2, width=width, val=val, err=err))
    
    def add_pressure(self,loc,val,err):
        '''! Add a constraint on total pressure

        @param loc (R,Z) position of constraint [m]
        @param val Pressure constraint [Pa]
        @param err Uncertainty in constraint [Pa]
        '''
        self._pressure_cons.append(Press_con(pt=loc,val=val,err=err))
    
    def reset_constraints(self):
        '''! Reset and remove all constraints'''
        self._Ip_con = None
        self._Dflux_con = None
        self._flux_loops = []
        self._mirnovs = []
        self._saddles = []
        self._pressure_cons = []
        self._coil_targets = None
        self._coil_wts = None
    
    def write_fit_in(self,filename='fit.in'):
        '''! Save constraints to a TokaMaker constraint file

        @param filename Path to constraint file
        '''
        constraints = self._flux_loops + self._mirnovs + self._pressure_cons
        if self._Ip_con is not None:
            constraints.append(self._Ip_con)
        if self._Dflux_con is not None:
            constraints.append(self._Dflux_con)
        ncons = len(constraints)
        with open(filename, 'w+') as fid:
            fid.write('{0:d}\n\n'.format(ncons))
            for con in constraints:
                con.write(fid)
    
    def read_fit_in(self,filename='fit.in'):
        '''! Read in constraints from a TokaMaker constraint file

        @param filename Path to constraint file
        '''
        self.reset_constraints()
        with open(filename, 'r') as fid:
            ncons = int(fid.readline())
            for _ in range(ncons):
                fid.readline()
                con_type = int(fid.readline())
                new_con_class = con_map[con_type]
                new_con = new_con_class()
                new_con.read(fid)
                if con_type == Ip_con_id:
                    self._Ip_con = new_con
                elif con_type == dFlux_con_id:
                    self._Dflux_con = new_con
                elif con_type == fluxLoop_con_id:
                    self._flux_loops.append(new_con)
                elif con_type == Mirnov_con_id:
                    self._mirnovs.append(new_con)
                elif con_type == saddle_con_id:
                    self._saddles.append(new_con)
                elif con_type == Press_con_id:
                    self._pressure_cons.append(new_con)
                else:
                    raise ValueError("Unknown constraint type")

    def reconstruct(self, vacuum=False):
        '''! Reconstruct G-S equation with specified fitting constraints, profiles, etc.
        
        @param vacuum Peform a vacuum reconstruction (no plasma current, pressure, etc.)
        '''
        self.write_fit_in()
        error_flag = c_int()
        if self._coil_targets is not None:
            self._gs_obj.set_coil_currents(self._coil_targets)
        if self._coil_wts is None:
            tokamaker_recon_run(c_bool(vacuum),self.settings,c_void_p(),ctypes.byref(error_flag))
        else:
            coil_wts_ptr = self._coil_wts.ctypes.data_as(c_void_p)
            tokamaker_recon_run(c_bool(vacuum),self.settings,coil_wts_ptr,ctypes.byref(error_flag))
        return error_flag.value