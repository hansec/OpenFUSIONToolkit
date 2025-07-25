# Auto-Generated on Fri Jan 24 11:50:18 2025
# using library build at /Users/jamiexia
# on machine JamBook-3.local
# settings: --nthread=2 --build_umfpack=1 --build_arpack=1 --oblas_threads

# Setup build and install paths
ROOT_PATH=$(pwd)
BUILD_DIR=$ROOT_PATH/build_release
INSTALL_DIR=$ROOT_PATH/install_release

# Create fresh build directory
rm -rf $BUILD_DIR
mkdir $BUILD_DIR && cd $BUILD_DIR

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR \
  -DOFT_BUILD_TESTS:BOOL=FALSE \
  -DOFT_BUILD_EXAMPLES:BOOL=TRUE \
  -DOFT_BUILD_PYTHON:BOOL=TRUE \
  -DOFT_BUILD_DOCS:BOOL=FALSE \
  -DOFT_USE_OpenMP:BOOL=TRUE \
  -DOFT_PACKAGE_BUILD:BOOL=FALSE \
  -DOFT_PACKAGE_NIGHTLY:BOOL=TRUE \
  -DOFT_COVERAGE:BOOL=FALSE \
  -DCMAKE_C_COMPILER:FILEPATH=gcc-14 \
  -DCMAKE_CXX_COMPILER:FILEPATH=g++-14 \
  -DCMAKE_Fortran_COMPILER:FILEPATH=gfortran-14 \
  -DCMAKE_Fortran_FLAGS:STRING="-fallow-argument-mismatch" \
  -DOFT_USE_MPI:BOOL=FALSE \
  -DOFT_METIS_ROOT:PATH=/Users/jamiexia/metis-5_1_0 \
  -DHDF5_ROOT:PATH=/Users/jamiexia/hdf5-1_14_5 \
  -DBLAS_ROOT:PATH=/Users/jamiexia/OpenBLAS-0_3_23 \
  -DLAPACK_ROOT:PATH=/Users/jamiexia/OpenBLAS-0_3_23 \
  -DBLA_VENDOR:STRING=OpenBLAS \
  -DOFT_ARPACK_ROOT:PATH=/Users/jamiexia/arpack-ng-3_5_0 \
  -DOFT_FoX_ROOT:PATH=/Users/jamiexia/fox-4_1_2 \
  -DOFT_UMFPACK_ROOT:PATH=/Users/jamiexia/UMFPACK-7_0_1 \
  /Users/jamiexia/OpenFUSIONToolkit/src
