################################################################################
#                                                                              #
#  MACHINE-SPECIFIC FUNCTIONS                                                  #
#                                                                              #
#    OPTIONS:                                                                  #
#      COMPILER   : PATH TO COMPILER EXECUTABLE                                #
#      GSL_DIR    : PATH TO GSL INSTALLATION                                   #
#      MPI_DIR    : PATH TO MPI INSTALLATION                                   #
#      HDF5_DIR   : PATH TO HDF5 INSTALLATION                                  #
#      EXECUTABLE : BINARY WRAPPER USED TO LAUNCH BHLIGHT                      #
#                                                                              #
#    MPI_DIR AND HDF5_DIR ARE NOT REQUIRED IF COMPILER HANDLES HEADERS AND     #
#    LIBRARIES FOR THESE DEPENDENCIES                                          #
#                                                                              #
################################################################################

import os

# module load phdf5
# module load gsl

def matches_host():
  host = os.uname()[1]
  return 'pfe' in host

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = '-O3 -fPIC -Wall -Werror -qopenmp -march=core-avx2'
  host['DEBUG_FLAGS']    = '-O0 -g -fPIC -Wall -Werror -qopenmp'
  host['HDF5_DIR']       = '/nasa/hdf5/1.8.18_mpt/'
  host['GSL_DIR']        = '/home4/ahankla/gsl/'
  host['EXECUTABLE']     = 'mpiexec'

  return host
