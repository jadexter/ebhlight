################################################################################
#                                                                              #
#  OPTICALLY THIN BREMSSTRAHLUNG COOLING                                       #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/');
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'gasbox'

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 6)
bhl.config.set_cparm('N2TOT', 6)
bhl.config.set_cparm('N3TOT', 6)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', False)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# ELECTRONS
bhl.config.set_cparm('ELECTRONS', True)
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', False)
bhl.config.set_cparm('BETA_HEAT', 3)
bhl.config.set_cparm('COULOMB', True)
bhl.config.set_cparm('TCOOL', 0)
bhl.config.set_cparm('INITELECTRONS', False)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'LINEAR')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')

# RADIATION
bhl.config.set_cparm('RADIATION', False)
bhl.config.set_cparm('COOLING', True)
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('NU_BINS_EMISS', 100)
bhl.config.set_cparm('NU_BINS_SPEC', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = 5120.e0)
bhl.config.set_rparm('dt', 'double', default = 5120.e-6)
bhl.config.set_rparm('tcool0', 'double', default=1.0)
bhl.config.set_rparm('tcoolOmega0', 'double', default=1.0)
bhl.config.set_rparm('q_constant', 'double', default=0.5)
bhl.config.set_rparm('L_unit', 'double', default = 1.)
bhl.config.set_rparm('M_unit', 'double', default = 1.)
# bhl.config.set_rparm('L_unit', 'double', default = 5.85532e7)
# bhl.config.set_rparm('M_unit', 'double', default = 2.01466e15)
bhl.config.set_rparm('tune_emiss', 'double', default = 1.e-5)
bhl.config.set_rparm('DTd', 'double', default = 5.)
bhl.config.set_rparm('DTl', 'double', default = 500.)
bhl.config.set_rparm('DTr', 'integer', default = 500000.)
bhl.config.set_rparm('Tp0', 'double', default = 1.e10)
bhl.config.set_rparm('Te0', 'double', default = 1.e8)
bhl.config.set_rparm('ne', 'double', default = 6.e19)
bhl.config.set_rparm('Tel_target', 'double', default=1e+9)
bhl.config.set_rparm('Tel_rslope', 'double', default=0.0)
# bhl.config.set_rparm('tp_over_te', 'double', default=3.0)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
