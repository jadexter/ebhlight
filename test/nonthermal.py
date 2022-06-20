################################################################################
#                                                                              #
# NONTHERMAL TESTING                                                           #
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
import glob
import numpy as np
import hdf5_to_dict as io
import util


TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'nonthermal'
AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
call(['python', 'build.py', '-dir', TMP_DIR])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
call(['./bhlight', '-p', 'param_template.dat'])
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
dump = io.load_dump(dfiles[-1], geom)
tscale = 1.e-2
gam = 1.4

rho_code = dump['RHO'][:,0,0]
P_code   = (gam-1.)*dump['UU'][:,0,0]/(tscale*tscale)
u_code   = dump['U1'][:,0,0]/(tscale)
#N1 = len(rho)
x_code = geom['x'][:,0,0]

# CLEAN UP
# util.safe_remove(TMP_DIR)

