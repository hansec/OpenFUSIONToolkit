#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''!Utility for manipulating and combining ThinCurr surface meshes.

This script operates on 3-node (linear) triangular surface meshes stored in the
Open FUSION Toolkit "native" HDF5 mesh format (as produced by @ref convert_cubit.py
"convert_cubit.py" and @ref convert_gmsh.py "convert_gmsh.py"). It supports two
top-level workflows via sub-commands:

Combine two (or more) existing mesh files into a single mesh:

    python ThinCurr_mesh_tool.py combine --in_files a.h5 b.h5 --out_file combined.h5

Modify a single mesh file by removing regions, applying a rigid transform, and/or
generating multiple shifted/rotated copies:

    # Remove regions 2 and 3
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --remove_regions 2 3

    # Shift the whole mesh by (0,0,0.5)
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --shift 0 0 0.5

    # Stretch the mesh by 2x along x (scale about the origin)
    python ThinCurr_mesh_tool.py modify --in_file a.h5 --scale 2 1 1

    # Build a 6-fold toroidal array (60 deg spacing about z) welding the seams:
    # the original plus 5 rotated copies at 60,120,180,240,300 deg
    python ThinCurr_mesh_tool.py modify --in_file segment.h5 --copies 5 \
        --rotate z 60 --weld_tol 1.0E-6

Only one of --shift, --rotate, or --scale may be given at a time. When --copies
is given, --shift or --rotate defines the per-copy increment (--scale is not
allowed); otherwise the transform is applied once to the whole mesh.

@note Only linear 3-node triangular meshes are supported. Meshes containing
high-order node information (``mesh/ho_info``) or non-triangular cells will be
rejected.

@authors Open FUSION Toolkit contributors
@date July 2026
'''
from __future__ import print_function
import argparse
import os
import numpy as np
import h5py


class ThinCurrMesh:
    '''!Container for a native-format 3-node triangular surface mesh.

    All connectivity is stored using 0-based indexing internally and converted
    to/from the 1-based convention used in the native HDF5 files on load/save.

    Attributes:
      r          Vertex coordinate list, shape [np,3] (float64); a 2D [np,2]
                 list supplied at construction is upgraded to 3D with Z=0
      lc         Cell (triangle) vertex list, shape [nc,3], 0-based (int32)
      reg        Per-cell region index, shape [nc], values in 1..nregions (int32)
      nodesets   List of node-index arrays (0-based); used for holes/closures
      sidesets   List of cell-index arrays (0-based)
      reg_attrs  Per-region attribute table, shape [nregions,nattr] or None
    '''

    def __init__(self, r, lc, reg, nodesets=None, sidesets=None, reg_attrs=None):
        '''!Construct a mesh from its component arrays

        @param r Vertex list [np,3]
        @param lc Triangle vertex list [nc,3], 0-based
        @param reg Per-cell region index [nc]
        @param nodesets List of 0-based node-index arrays (optional)
        @param sidesets List of 0-based cell-index arrays (optional)
        @param reg_attrs Per-region attribute table [nregions,nattr] (optional)
        '''
        self.r = np.ascontiguousarray(r, dtype=np.float64)
        self.lc = np.ascontiguousarray(lc, dtype=np.int32)
        self.reg = np.ascontiguousarray(reg, dtype=np.int32).reshape(-1)
        if self.lc.ndim != 2 or self.lc.shape[1] != 3:
            raise ValueError("Cell list must have shape [nc,3] (3-node triangles)")
        if self.r.ndim != 2 or self.r.shape[1] not in (2, 3):
            raise ValueError("Vertex list must have shape [np,2] or [np,3]")
        if self.r.shape[1] == 2:  # Upgrade 2D point list to 3D (Z=0)
            self.r = np.hstack((self.r, np.zeros((self.r.shape[0], 1))))
        if self.reg.shape[0] != self.lc.shape[0]:
            raise ValueError("Region list length must match number of cells")
        self.nodesets = [np.asarray(ns, dtype=np.int32).reshape(-1) for ns in (nodesets or [])]
        self.sidesets = [np.asarray(ss, dtype=np.int32).reshape(-1) for ss in (sidesets or [])]
        if reg_attrs is None:
            self.reg_attrs = None
        else:
            self.reg_attrs = np.ascontiguousarray(reg_attrs, dtype=np.float64)
            if self.reg_attrs.ndim != 2:
                raise ValueError("Region attributes must be a 2D table [nregions,nattr]")

    @property
    def np(self):
        '''!Number of vertices'''
        return self.r.shape[0]

    @property
    def nc(self):
        '''!Number of cells (triangles)'''
        return self.lc.shape[0]

    @property
    def nregions(self):
        '''!Number of regions (maximum region index)'''
        if self.nc == 0:
            return 0
        return int(self.reg.max())

    def copy(self):
        '''!Return a deep copy of the mesh'''
        return ThinCurrMesh(
            self.r.copy(), self.lc.copy(), self.reg.copy(),
            [ns.copy() for ns in self.nodesets],
            [ss.copy() for ss in self.sidesets],
            None if self.reg_attrs is None else self.reg_attrs.copy())

    # ------------------------------------------------------------------ I/O
    @classmethod
    def load(cls, filename):
        '''!Load a mesh from a native-format HDF5 file

        @param filename Path to input mesh file
        @result New @ref ThinCurrMesh instance
        '''
        print()
        print("Reading mesh: {0}".format(filename))
        with h5py.File(filename, 'r') as fid:
            if 'mesh/ho_info' in fid:
                raise ValueError("High-order meshes are not supported (found 'mesh/ho_info')")
            r = np.asarray(fid['mesh/R'])
            lc = np.asarray(fid['mesh/LC']) - 1  # convert to 0-based
            reg = np.asarray(fid['mesh/REG'])
            if lc.shape[1] != 3:
                raise ValueError("Only 3-node triangular meshes are supported "
                                 "(found {0} nodes/cell)".format(lc.shape[1]))
            if r.ndim == 2 and r.shape[1] == 2:
                print("  Note: input mesh is 2D; padding point list to 3D (Z=0)")
            nodesets = []
            if 'mesh/NUM_NODESETS' in fid:
                for j in range(int(fid['mesh/NUM_NODESETS'][0])):
                    nodesets.append(np.asarray(fid['mesh/NODESET{0:04d}'.format(j + 1)]) - 1)
            sidesets = []
            if 'mesh/NUM_SIDESETS' in fid:
                for j in range(int(fid['mesh/NUM_SIDESETS'][0])):
                    sidesets.append(np.asarray(fid['mesh/SIDESET{0:04d}'.format(j + 1)]) - 1)
            reg_attrs = None
            if 'mesh/reg_attr/NUM_ATTR' in fid:
                nattr = int(fid['mesh/reg_attr/NUM_ATTR'][0])
                cols = [np.asarray(fid['mesh/reg_attr/ATTR{0:04d}'.format(k + 1)]) for k in range(nattr)]
                reg_attrs = np.vstack(cols).T if nattr > 0 else None
        mesh = cls(r, lc, reg, nodesets, sidesets, reg_attrs)
        mesh._check_attrs()
        mesh.print_info()
        return mesh

    def save(self, filename):
        '''!Save the mesh to a native-format HDF5 file

        @param filename Path to output mesh file
        '''
        self._check_attrs()
        print()
        print("Saving mesh: {0}".format(filename))
        with h5py.File(filename, 'w') as fid:
            fid.create_dataset('mesh/R', data=self.r, dtype='f8')
            fid.create_dataset('mesh/LC', data=self.lc + 1, dtype='i4')  # 1-based
            fid.create_dataset('mesh/REG', data=self.reg, dtype='i4')
            if self.reg_attrs is not None and self.reg_attrs.shape[1] > 0:
                nattr = self.reg_attrs.shape[1]
                fid.create_dataset('mesh/reg_attr/NUM_ATTR', data=[nattr], dtype='i4')
                for k in range(nattr):
                    fid.create_dataset('mesh/reg_attr/ATTR{0:04d}'.format(k + 1),
                                       data=self.reg_attrs[:, k], dtype='f8')
            if len(self.nodesets) > 0:
                fid.create_dataset('mesh/NUM_NODESETS', data=[len(self.nodesets)], dtype='i4')
                for j, ns in enumerate(self.nodesets):
                    fid.create_dataset('mesh/NODESET{0:04d}'.format(j + 1), data=ns + 1, dtype='i4')
            if len(self.sidesets) > 0:
                fid.create_dataset('mesh/NUM_SIDESETS', data=[len(self.sidesets)], dtype='i4')
                for j, ss in enumerate(self.sidesets):
                    fid.create_dataset('mesh/SIDESET{0:04d}'.format(j + 1), data=ss + 1, dtype='i4')
        self.print_info()

    def print_info(self):
        '''!Print a short summary of the mesh contents'''
        print("  # of points   = {0}".format(self.np))
        print("  # of cells    = {0}".format(self.nc))
        print("  # of regions  = {0}".format(self.nregions))
        print("  # of nodesets = {0}".format(len(self.nodesets)))
        print("  # of sidesets = {0}".format(len(self.sidesets)))
        if self.reg_attrs is not None:
            print("  # of attrs    = {0}".format(self.reg_attrs.shape[1]))

    def _check_attrs(self):
        '''!Validate that the attribute table (if present) matches the region count'''
        if self.reg_attrs is not None and self.reg_attrs.shape[0] != self.nregions:
            raise ValueError("Region attribute table has {0} rows but mesh has {1} regions".format(
                self.reg_attrs.shape[0], self.nregions))

    # ------------------------------------------------------------- editing
    def transform(self, shift=None, rotate=None, scale=None, center=None):
        '''!Apply a transform to the vertex coordinates (in place)

        Operations, if supplied, are applied in the order rotate, scale, shift:
        ``r' = S R (r - center) + center + shift`` (scaling is about the origin).

        @param shift Translation vector [3] (optional)
        @param rotate Tuple ``(axis, angle_deg)`` where axis is 'x', 'y', or 'z' (optional)
        @param scale Per-axis scale factors [3] applied about the origin (optional)
        @param center Center of rotation [3] (default: origin)
        '''
        if rotate is not None:
            axis, angle_deg = rotate
            rmat = rotation_matrix(axis, angle_deg)
            c = np.zeros(3) if center is None else np.asarray(center, dtype=np.float64)
            self.r = (self.r - c) @ rmat.T + c
        if scale is not None:
            self.r = self.r * np.asarray(scale, dtype=np.float64)
        if shift is not None:
            self.r = self.r + np.asarray(shift, dtype=np.float64)
        return self

    def remove_regions(self, regions):
        '''!Remove all cells belonging to the specified regions (in place)

        Surviving regions are renumbered contiguously (1..N), unreferenced
        vertices are dropped, and nodesets/sidesets are remapped, dropping any
        entries that reference removed nodes/cells.

        @param regions Iterable of region indices to remove
        '''
        remove = set(int(x) for x in regions)
        missing = remove - set(int(x) for x in np.unique(self.reg))
        if missing:
            print("  Warning: region(s) {0} not present in mesh".format(sorted(missing)))
        keep_cell = np.array([r not in remove for r in self.reg], dtype=bool)
        n_removed = int(np.sum(~keep_cell))
        print("  Removing {0} cell(s) in region(s) {1}".format(n_removed, sorted(remove)))
        # Renumber surviving regions to a contiguous 1..N range
        surviving = np.unique(self.reg[keep_cell]) if np.any(keep_cell) else np.array([], dtype=np.int32)
        remap = {int(old): new + 1 for new, old in enumerate(surviving)}
        if self.reg_attrs is not None and surviving.shape[0] > 0:
            self.reg_attrs = self.reg_attrs[surviving - 1, :].copy()
        elif self.reg_attrs is not None:
            self.reg_attrs = self.reg_attrs[0:0, :].copy()
        self.lc = self.lc[keep_cell, :]
        self.reg = np.array([remap[int(r)] for r in self.reg[keep_cell]], dtype=np.int32)
        # Remap sidesets (cell indices)
        cell_new = -np.ones((keep_cell.shape[0],), dtype=np.int64)
        cell_new[keep_cell] = np.arange(int(np.sum(keep_cell)))
        self.sidesets = _remap_index_sets(self.sidesets, cell_new, "sideset")
        self._reindex_vertices()
        return self

    def weld(self, tol):
        '''!Merge coincident vertices within a distance tolerance (in place)

        Vertices whose rounded coordinates (to the nearest `tol`) coincide are
        merged to a single node. Degenerate triangles produced by welding (two
        or more identical corners) are removed. Useful for stitching seams
        between replicated or combined mesh copies.

        @param tol Merge tolerance (grid spacing) applied to each coordinate
        '''
        if tol <= 0.0:
            raise ValueError("Weld tolerance must be positive")
        keys = np.round(self.r / tol).astype(np.int64)
        _, first_idx, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
        inverse = inverse.reshape(-1)
        n_merged = self.np - first_idx.shape[0]
        if n_merged == 0:
            print("  Weld: no coincident vertices found (tol={0})".format(tol))
            return self
        # New vertex ordering follows first occurrence to keep results stable
        order = np.argsort(first_idx)
        new_of_unique = np.empty((first_idx.shape[0],), dtype=np.int64)
        new_of_unique[order] = np.arange(first_idx.shape[0])
        old_to_new = new_of_unique[inverse]
        self.r = self.r[first_idx[order], :]
        self.lc = old_to_new[self.lc].astype(np.int32)
        self.nodesets = [np.unique(old_to_new[ns]).astype(np.int32) for ns in self.nodesets]
        # Drop degenerate triangles created by welding
        good = (self.lc[:, 0] != self.lc[:, 1]) & (self.lc[:, 1] != self.lc[:, 2]) & (self.lc[:, 0] != self.lc[:, 2])
        n_degen = int(np.sum(~good))
        if n_degen > 0:
            cell_new = -np.ones((good.shape[0],), dtype=np.int64)
            cell_new[good] = np.arange(int(np.sum(good)))
            self.lc = self.lc[good, :]
            self.reg = self.reg[good]
            self.sidesets = _remap_index_sets(self.sidesets, cell_new, "sideset")
        print("  Weld: merged {0} vertices, removed {1} degenerate cell(s) (tol={2})".format(
            n_merged, n_degen, tol))
        self.nodesets = [ns for ns in self.nodesets if ns.shape[0] > 0]
        return self

    def _reindex_vertices(self):
        '''!Drop unreferenced vertices and remap connectivity/nodesets (in place)'''
        used = np.zeros((self.np,), dtype=bool)
        if self.nc > 0:
            used[self.lc.reshape(-1)] = True
        new_of_old = -np.ones((self.np,), dtype=np.int64)
        new_of_old[used] = np.arange(int(np.sum(used)))
        self.r = self.r[used, :]
        if self.nc > 0:
            self.lc = new_of_old[self.lc].astype(np.int32)
        self.nodesets = _remap_index_sets(self.nodesets, new_of_old, "nodeset")

    def append(self, other, distinct_regions=True):
        '''!Append another mesh onto this one (in place)

        @param other The @ref ThinCurrMesh to append
        @param distinct_regions If True, offset the appended mesh's region indices
            so its regions remain distinct from this mesh's regions
        @result self
        '''
        np_offset = self.np
        nc_offset = self.nc
        reg_offset = self.nregions if distinct_regions else 0
        self.r = np.vstack((self.r, other.r)) if self.np > 0 else other.r.copy()
        self.lc = np.vstack((self.lc, other.lc + np_offset)) if nc_offset > 0 else (other.lc + np_offset)
        self.reg = np.concatenate((self.reg, other.reg + reg_offset))
        for ns in other.nodesets:
            self.nodesets.append(ns + np_offset)
        for ss in other.sidesets:
            self.sidesets.append(ss + nc_offset)
        # Attributes: carry through only if compatible and regions kept distinct
        self.reg_attrs = _combine_attrs(self.reg_attrs, other.reg_attrs,
                                        reg_offset > 0 or self.nregions == other.nregions)
        return self


def rotation_matrix(axis, angle_deg):
    '''!Build a rotation matrix about a principal Cartesian axis

    @param axis Rotation axis: 'x', 'y', or 'z' (case-insensitive)
    @param angle_deg Rotation angle in degrees (right-handed about `axis`)
    @result 3x3 rotation matrix
    '''
    a = str(axis).lower()
    t = np.deg2rad(float(angle_deg))
    c, s = np.cos(t), np.sin(t)
    if a == 'x':
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=np.float64)
    if a == 'y':
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=np.float64)
    if a == 'z':
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=np.float64)
    raise ValueError("Rotation axis must be one of 'x', 'y', or 'z' (got '{0}')".format(axis))


def _remap_index_sets(sets, old_to_new, label):
    '''!Remap and filter a list of index arrays using an old->new lookup

    Entries mapped to a negative value (removed) are dropped. Empty sets that
    result are discarded with a warning.

    @param sets List of index arrays (0-based)
    @param old_to_new Lookup array mapping old index -> new index (<0 = removed)
    @param label Name used in warning messages
    @result Filtered/remapped list of index arrays
    '''
    out = []
    for i, s in enumerate(sets):
        mapped = old_to_new[s]
        keep = mapped >= 0
        if not np.all(keep):
            print("  Warning: dropping {0} entrie(s) from {1} {2}".format(
                int(np.sum(~keep)), label, i + 1))
        mapped = mapped[keep].astype(np.int32)
        if mapped.shape[0] > 0:
            out.append(mapped)
        else:
            print("  Warning: {0} {1} became empty and was dropped".format(label, i + 1))
    return out


def _combine_attrs(a, b, compatible):
    '''!Stack two per-region attribute tables along the region axis

    @param a First attribute table or None
    @param b Second attribute table or None
    @param compatible Whether the two meshes' regions are being kept distinct
    @result Combined attribute table or None (with a warning when dropped)
    '''
    if a is None and b is None:
        return None
    if a is None or b is None:
        print("  Warning: region attributes present on only one mesh; dropping attributes")
        return None
    if a.shape[1] != b.shape[1] or not compatible:
        print("  Warning: incompatible region attributes; dropping attributes")
        return None
    return np.vstack((a, b))


def combine_meshes(filenames, distinct_regions=True, weld_tol=None):
    '''!Combine multiple mesh files into a single mesh

    @param filenames List of input file paths
    @param distinct_regions Keep each input's regions distinct by offsetting IDs
    @param weld_tol Optional tolerance for welding coincident vertices at seams
    @result Combined @ref ThinCurrMesh
    '''
    meshes = [ThinCurrMesh.load(fn) for fn in filenames]
    combined = meshes[0].copy()
    for mesh in meshes[1:]:
        combined.append(mesh, distinct_regions=distinct_regions)
    if weld_tol is not None:
        combined.weld(weld_tol)
    return combined


def replicate_mesh(mesh, ncopies, shift=None, rotate=None, center=None,
                   distinct_regions=True, weld_tol=None):
    '''!Generate an output mesh consisting of the original plus transformed copies

    The untransformed original is always kept. Copy ``k`` (k = 1 .. ncopies) is
    the input mesh with the incremental transform applied ``k`` times: shifted by
    ``k*shift`` and/or rotated by ``k*angle`` about `center`. The output therefore
    contains ``ncopies + 1`` instances.

    @param mesh Input @ref ThinCurrMesh
    @param ncopies Number of transformed copies to add (excludes the original)
    @param shift Per-copy incremental translation [3] (optional)
    @param rotate Per-copy incremental rotation ``(axis, angle_deg)`` (optional)
    @param center Center of rotation [3] (default: origin)
    @param distinct_regions If True, offset region IDs so every copy is distinct
    @param weld_tol Optional tolerance for welding coincident vertices at seams
    @result Replicated @ref ThinCurrMesh
    '''
    if ncopies < 1:
        raise ValueError("Number of copies must be >= 1")
    if shift is None and rotate is None:
        raise ValueError("Replication requires a --shift or --rotate increment")
    out = None
    for i in range(ncopies + 1):  # i=0 is the untransformed original
        piece = mesh.copy()
        step_shift = None if shift is None else np.asarray(shift, dtype=np.float64) * i
        step_rotate = None if rotate is None else (rotate[0], float(rotate[1]) * i)
        piece.transform(shift=step_shift, rotate=step_rotate, center=center)
        if out is None:
            out = piece
        else:
            out.append(piece, distinct_regions=distinct_regions)
    if weld_tol is not None:
        out.weld(weld_tol)
    return out


def _add_common_out(parser):
    parser.add_argument("--out_file", type=str, default=None, help="Output mesh file")


def build_parser():
    '''!Construct the command-line argument parser'''
    parser = argparse.ArgumentParser(
        description="Manipulate and combine native-format 3-node triangular ThinCurr meshes")
    sub = parser.add_subparsers(dest="command", required=True)

    # --- combine ---
    p_comb = sub.add_parser("combine", help="Combine two or more mesh files")
    p_comb.add_argument("--in_files", type=str, nargs='+', required=True,
                        help="Input mesh files (two or more)")
    p_comb.add_argument("--merge_regions", action="store_true", default=False,
                        help="Merge identical region IDs across inputs instead of "
                             "keeping them distinct (default: keep distinct)")
    p_comb.add_argument("--weld_tol", type=float, default=None,
                        help="Weld coincident vertices within this tolerance")
    _add_common_out(p_comb)

    # --- modify ---
    p_mod = sub.add_parser("modify", help="Remove regions, transform, and/or replicate a mesh")
    p_mod.add_argument("--in_file", type=str, required=True, help="Input mesh file")
    p_mod.add_argument("--remove_regions", type=int, nargs='+', default=None,
                       help="Region indices to remove")
    xform = p_mod.add_mutually_exclusive_group()
    xform.add_argument("--shift", type=float, nargs=3, default=None,
                       metavar=('X', 'Y', 'Z'),
                       help="Translation [X Y Z]. Applied once to the whole mesh, or as "
                            "the per-copy increment when --copies is given")
    xform.add_argument("--rotate", nargs=2, default=None, metavar=('AXIS', 'ANGLE'),
                       help="Rotation of ANGLE degrees about AXIS (x|y|z). Applied once to "
                            "the whole mesh, or as the per-copy increment when --copies is given")
    xform.add_argument("--scale", type=float, nargs=3, default=None,
                       metavar=('SX', 'SY', 'SZ'),
                       help="Scale factors along X, Y, Z about the origin. Cannot be "
                            "combined with --copies")
    p_mod.add_argument("--rotate_center", type=float, nargs=3, default=None,
                       metavar=('X', 'Y', 'Z'), help="Center of rotation (default: origin)")
    p_mod.add_argument("--copies", type=int, default=None,
                       help="Add this many transformed copies, not counting the original "
                            "(which is always kept), applying --shift or --rotate "
                            "incrementally to each (requires --shift or --rotate; not "
                            "valid with --scale)")
    p_mod.add_argument("--persist_regions", action="store_true", default=False,
                       help="Keep the original region IDs on every copy instead of "
                            "offsetting them so each copy is distinct (default: offset "
                            "so each copy gets its own distinct regions)")
    p_mod.add_argument("--weld_tol", type=float, default=None,
                       help="Weld coincident vertices within this tolerance")
    _add_common_out(p_mod)
    return parser


def _parse_rotate(val):
    '''!Validate and normalize a (axis, angle) argument pair'''
    if val is None:
        return None
    axis, angle = val
    rotation_matrix(axis, angle)  # validates axis and angle
    return (str(axis).lower(), float(angle))


def main(argv=None):
    '''!Command-line entry point'''
    parser = build_parser()
    options = parser.parse_args(argv)

    if options.command == "combine":
        if len(options.in_files) < 2:
            parser.error("combine requires at least two --in_files")
        out_file = options.out_file
        if out_file is None:
            out_file = os.path.splitext(options.in_files[0])[0] + "-combined.h5"
        result = combine_meshes(options.in_files,
                                distinct_regions=not options.merge_regions,
                                weld_tol=options.weld_tol)
        result.save(out_file)

    elif options.command == "modify":
        try:
            rotate = _parse_rotate(options.rotate)
        except ValueError as err:
            parser.error(str(err))
        # --shift/--rotate/--scale are mutually exclusive via argparse; only the
        # additional constraints tied to --copies need checking here.
        if options.copies is not None:
            if options.scale is not None:
                parser.error("--scale cannot be combined with --copies")
            if options.shift is None and rotate is None:
                parser.error("--copies requires --shift or --rotate")
        out_file = options.out_file
        if out_file is None:
            out_file = os.path.splitext(options.in_file)[0] + "-modified.h5"

        mesh = ThinCurrMesh.load(options.in_file)
        if options.remove_regions is not None:
            print()
            print("Removing regions...")
            mesh.remove_regions(options.remove_regions)
        if options.copies is not None:
            print()
            print("Replicating mesh...")
            mesh = replicate_mesh(mesh, options.copies,
                                  shift=options.shift, rotate=rotate,
                                  center=options.rotate_center,
                                  distinct_regions=not options.persist_regions,
                                  weld_tol=options.weld_tol)
        else:
            if options.shift is not None or rotate is not None or options.scale is not None:
                print()
                print("Applying transform...")
                mesh.transform(shift=options.shift, rotate=rotate,
                               scale=options.scale, center=options.rotate_center)
            if options.weld_tol is not None:
                mesh.weld(options.weld_tol)
        mesh.save(out_file)


if __name__ == "__main__":
    main()
