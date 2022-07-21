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
  return 'shas0430.rc.int.colorado.edu' in host

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = '-O3 -fPIC -Wall -Werror -fopenmp'
  host['DEBUG_FLAGS']    = '-O0 -g -fPIC -Wall -Werror -fopenmp'
  host['HDF5_DIR']       = '/curc/sw/hdf5/1.12.1/gcc/10.2.0'
  host['GSL_DIR']        = '/curc/sw/gsl/2.6/gcc/10.2.0'
  host['EXECUTABLE']     = 'mpirun'

  return host

