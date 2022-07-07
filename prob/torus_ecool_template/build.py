import sys; sys.path.append('../../script');
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'torus'

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 320)
bhl.config.set_cparm('N2TOT', 256)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', 10) # Match these to job scripts
bhl.config.set_cparm('N2CPU', 8)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MMKS')

# ELECTRONS
bhl.config.set_cparm('ELECTRONS', True)
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', False) # not actually used anywhere
bhl.config.set_cparm('BETA_HEAT', 3) # Sets electron heating prescription.
bhl.config.set_cparm('COULOMB', True) # not actually used anywhere (i.e. can't turn off)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

# COOLING
bhl.config.set_cparm('COOLING', True) # electron-only cooling function
bhl.config.set_cparm('TCOOL', 1) # Sets electron cooling time. Only relevant if COOLING. Option 1 sets tcool = 1/Omega. TCOOL 0 sets tcool to a constant (see tcool0 below)
bhl.config.set_cparm('INITELECTRONS', True)

# RADIATION
bhl.config.set_cparm('RADIATION', False)
# The following parameters only come into effect if RADIATION
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', True)
bhl.config.set_cparm('NU_BINS_EMISS', 100)
bhl.config.set_cparm('NU_BINS_SPEC', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('SYNCHROTRON', True)
# Radiation boundary conditions
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_CAMERA')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###
# Electron-cooling function parameters
bhl.config.set_rparm('Tel_target', 'double', default=1e+9) # Target electron temperature
bhl.config.set_rparm('Tel_rslope', 'double', default=0.0) # Power-law index for radial temp dependence
bhl.config.set_rparm('tcool0', 'double', default=1.0) # Only used if TCOOL is 0; fixes tcool to set value instead of 1/Omega
bhl.config.set_rparm('tcoolOmega0', 'double', default=1.0) # Only used if TCOOL is 1; fixes tcool=1/(Omega*tcoolOmega0)
bhl.config.set_rparm('q_constant', 'double', default=0.5) # Only used if TCOOL is 1; Noble+ uses 0.5, Nico uses 1.
# Standard run-time parameters
bhl.config.set_rparm('tf', 'double', default = 2501.0) # Simulation end time
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('Rout', 'double', default = 1000.) # Outer extent of simulation domain
bhl.config.set_rparm('Rout_rad', 'double', default = 40.) # only if RADIATION
bhl.config.set_rparm('gam', 'double', default = 13./9.) # Total gas
bhl.config.set_rparm('DTd', 'double', default = 5.) # dump frequency
bhl.config.set_rparm('DTl', 'double', default = 5.e-1) # log frequency
bhl.config.set_rparm('DTr', 'double', default = 1000.) # restart frequency
bhl.config.set_rparm('DNr', 'integer', default = 1000)
# Other problem-specific parameters
bhl.config.set_rparm('rin', 'double', default = 20)     # Inner edge of torus
bhl.config.set_rparm('rmax', 'double', default = 41)    # Torus pressure maximum
bhl.config.set_rparm('a', 'double', default = 0.9375) # black hole spin
bhl.config.set_rparm('mbh', 'double', default = 1.e8) # black hole mass (in M_unit)
bhl.config.set_rparm('M_unit', 'double', default = 8.e23) # mass unit
bhl.config.set_rparm('u_jitter', 'double', default = 0.04) # fluid perturbations
bhl.config.set_rparm('MAD', 'int', default = 1)
bhl.config.set_rparm('BHflux', 'double', default = 0.)
bhl.config.set_rparm('beta', 'double', default = 100.)
bhl.config.set_rparm('tp_over_te', 'double', default=3) # only if not ELECTRONS
bhl.config.set_rparm('cour', 'double', default=0.9) # courant condition
# Below parameters only matter if RADIATION
bhl.config.set_rparm('tune_emiss', 'double', 1.0)
bhl.config.set_rparm('tune_scatt', 'double', 0.1)
bhl.config.set_rparm('t0_tune_emiss', 'double', 500)
bhl.config.set_rparm('t0_tune_scatt', 'double', 500)
bhl.config.set_rparm('numin_emiss', 'double', default=1.e10)
bhl.config.set_rparm('numax_emiss', 'double', default=1.e16)
bhl.config.set_rparm('numin_spec', 'double', default=1.e10)
bhl.config.set_rparm('numax_spec', 'double', default=1.e25)
bhl.config.set_rparm('nph_per_proc', 'double', default=1.e5)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
