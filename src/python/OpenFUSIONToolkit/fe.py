#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Access to finite element functionality for Open FUSION Toolkit (OFT) Python interfaces

@authors Chris Hansen
@date September 2026
@ingroup doxy_oft_python
'''
import ctypes
import numpy
from ._interface import *


class FE_field_interpolator():
    '''! Base interpolation class for FE fields

    Wraps a Fortran field evaluation object, created by @ref Lagrange_surface_field_interpolator "a generic constructor"
    or a tool-specific constructor (eg. @ref OpenFUSIONToolkit.Marklin.Marklin.get_binterp "Marklin.get_binterp()"),
    and evaluates it using the generic FE evaluation interface. The Fortran object is destroyed
    when this object is deleted.
    '''
    def __init__(self,int_obj,dim,fbary_tol=1.E-8,dim_eval=None):
        '''! Initialize interpolation object

        @param int_obj Address of Fortran field evaluation object
        @param dim Dimension of field returned by @ref eval
        @param fbary_tol Tolerance for physical to logical mapping
        @param dim_eval Dimension of field returned by Fortran evaluation (default: `dim`)
        '''
        ## Address of Fortran field evaluation object
        self._int_obj = int_obj
        ## Dimension of field
        self.dim = dim
        ## Dimension of field returned by Fortran evaluation (may be larger than `dim`)
        self._dim_eval = dim if dim_eval is None else dim_eval
        ## Tolerance for physical to logical mapping
        self.fbary_tol = fbary_tol
        ## Cell containing the last point evaluated (0: not found; negative: outside of `fbary_tol`),
        ## which is used as the initial guess for the next search
        self.cell = c_int(0)

    def __del__(self):
        '''Destroy underlying interpolation object'''
        int_obj = getattr(self,'_int_obj',None)
        if int_obj:
            oft_destroy_field_eval(int_obj)
            self._int_obj = c_void_p()

    def _eval(self,pts):
        '''! Evaluate field at 3D locations

        @param pts Locations for evaluation [n,3]
        @result Field at evaluation points [n,self._dim_eval] (NaN for points outside the mesh)
        '''
        pts = numpy.ascontiguousarray(pts, dtype=numpy.float64)
        npts = pts.shape[0]
        vals_out = numpy.zeros((npts,self._dim_eval), dtype=numpy.float64)
        oft_apply_field_eval(self._int_obj,pts,npts,self.fbary_tol,self._dim_eval,ctypes.byref(self.cell),vals_out)
        return vals_out

    def eval(self,pts):
        '''! Evaluate field at given locations

        @param pts Location(s) for evaluation [3] or [n,3]
        @result Field at evaluation point(s) [self.dim] or [n,self.dim] (NaN for points outside the mesh)
        '''
        pts = numpy.asarray(pts, dtype=numpy.float64)
        in_1d = (pts.ndim == 1)
        if in_1d:
            pts = pts.reshape((1,pts.shape[0]))
        if pts.shape[1] != 3:
            raise ValueError('Input points must be 3D coordinates')
        vals_out = self._eval(pts)
        if in_1d:
            return vals_out[0,:self.dim].copy()
        else:
            return vals_out[:,:self.dim].copy()


class Lagrange_surface_field_interpolator(FE_field_interpolator):
    '''! Interpolation class for surface Lagrange FE fields'''
    def __init__(self,oft_env,fe_obj,values,int_type,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param oft_env OFT runtime environment object (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        @param fe_obj Address of Fortran FE object
        @param values FE weight vector
        @param int_type Interpolation type (1: value; 2: gradient)
        @param fbary_tol Tolerance for physical to logical mapping
        '''
        if int_type == 1:
            dim = 1
            imode = 211
        elif int_type == 2:
            dim = 3
            imode = 212
        else:
            raise ValueError('Invalid interpolation type, must be 1 or 2')
        self._oft_env = oft_env
        self._fe_obj = fe_obj
        #
        values = numpy.ascontiguousarray(values, dtype=numpy.float64)
        int_obj = c_void_p()
        error_string = self._oft_env.get_c_errorbuff()
        oft_get_field_eval(fe_obj,values,imode,ctypes.byref(int_obj),error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        super().__init__(int_obj,dim,fbary_tol=fbary_tol)


class Lagrange_2D_field_interpolator(Lagrange_surface_field_interpolator):
    '''! Interpolation class for 2D Lagrange FE fields'''
    def __init__(self,oft_env,fe_obj,values,int_type,fbary_tol=1.E-8):
        '''! Initialize interpolation object

        @param oft_env OFT runtime environment object (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        @param fe_obj Address of Fortran FE object
        @param values FE weight vector
        @param int_type Interpolation type (1: value; 2: gradient)
        @param fbary_tol Tolerance for physical to logical mapping
        '''
        super().__init__(oft_env,fe_obj,values,int_type,fbary_tol)
        if self.dim == 3: # Drop z-component for 2D interpolation
            self.dim = 2

    def eval(self,pts):
        '''! Evaluate field at given locations

        @param pts Location(s) for evaluation [2] or [n,2]
        @result Field at evaluation point(s) [self.dim] or [n,self.dim] (NaN for points outside the mesh)
        '''
        pts = numpy.asarray(pts, dtype=numpy.float64)
        in_1d = (pts.ndim == 1)
        if in_1d:
            pts = pts.reshape((1,pts.shape[0]))
        if pts.shape[1] != 2:
            raise ValueError('Input points must be 2D coordinates')
        pts_eval = numpy.zeros((pts.shape[0],3), dtype=numpy.float64)
        pts_eval[:,:2] = pts
        vals_out = self._eval(pts_eval)
        if in_1d:
            return vals_out[0,:self.dim].copy()
        else:
            return vals_out[:,:self.dim].copy()
