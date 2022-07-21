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

import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import root


TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'nonthermal_basic'
AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

##### GENERAL TESTING TOOLS #####

def syngam0(gamma,t):
    return 1./((1./gamma)-bs*t)

def syngam1(gammax,t):
    return 1./((1./gammax)+bs*t)

def syngam2(gammin,t):
    return 1./((1./gammin)+bs*t)

def synTau(gammaprime, gamma):
    return (bs*gammaprime)**-1 - (bs*gamma)**-1

def adiabgam0(gamma,t):
    return np.sqrt(((gamma**2)-1)*np.exp(2*bad*t)+1)

def adiabgam1(gammin,t):
    return np.sqrt(((gammin**2)-1)*np.exp(-2*bad*t)+1)

def adiabgam2(gammax,t):
    return np.sqrt(((gammax**2)-1)*np.exp(-2*bad*t)+1)

def adiabTau(gammaprime,gamma):
    return (1./(2*bad))*np.log(((gammaprime**2)-1)/((gamma**2)-1))

def Q(gamma):
    return Q0*(gamma**-p)

def step(x):
    return 1 * (x > 0)

def integrand(x,t,gamma):
    return Q(x)

def N(t,gammalist,radtype):
    if radtype=='syn':
        gamma1 = syngam1(ginjmax,t)
        gamma2 = syngam2(ginjmin,t)
        test1 = (ginjmin < gamma1)
        test2 = (gamma1 <= ginjmin)

        numgam = np.zeros(np.size(gammalist))

        i = 0
        for gamma in gammalist:
            if (gamma2<gamma)&(gamma<=ginjmax):
                if ginjmin < gamma1:
                    if (gamma<=ginjmin):
                        numgam[i],foo = quad(integrand,ginjmin,syngam0(gamma,t),args=(t,gamma))
                    elif (ginjmin<gamma)&(gamma<=gamma1):
                        numgam[i],foo = quad(integrand,gamma,syngam0(gamma,t),args=(t,gamma))
                    else:
                        numgam[i],foo = quad(integrand,gamma,ginjmax,args=(t,gamma))
                else:
                    if gamma <= gamma1:
                        numgam[i],foo = quad(integrand,ginjmin,syngam0(gamma,t),args=(t,gamma))
                    elif (gamma1<gamma)&(gamma<=ginjmin):
                        numgam[i],foo = quad(integrand,ginjmin,ginjmax,args=(t,gamma))
                    else:
                        numgam[i],foo = quad(integrand,gamma,ginjmax,args=(t,gamma))


            numgam[i] = (1/(bs*gamma**2))*numgam[i]
            i+=1
            
            
    elif radtype=='adiab':
        gamma1 = adiabgam1(ginjmin,t)
        gamma2 = adiabgam2(ginjmax,t)

        numgam = np.zeros(np.size(gammalist))

        i = 0
        for gamma in gammalist:
            if (ginjmin<gamma)&(gamma<=gamma2):
                if gamma1 < ginjmax:
                    if (ginjmax<gamma):
                        numgam[i],foo = quad(integrand,adiabgam0(gamma,t),ginjmax,args=(t,gamma))
                    elif (gamma1<gamma)&(gamma<=ginjmax):
                        numgam[i],foo = quad(integrand,adiabgam0(gamma,t),gamma,args=(t,gamma))
                    else:
                        numgam[i],foo = quad(integrand,ginjmin,gamma,args=(t,gamma))
                else:
                    if gamma1 < gamma:
                        numgam[i],foo = quad(integrand,adiabgam0(gamma,t),ginjmax,args=(t,gamma))
                    elif (ginjmax<gamma)&(gamma<=gamma1):
                        numgam[i],foo = quad(integrand,ginjmin,ginjmax,args=(t,gamma))
                    else:
                        numgam[i],foo = quad(integrand,ginjmin,gamma,args=(t,gamma))
            
            # neg sign is to flip integral (ie. lower is the smallest gamma that could have heated)
            numgam[i] = -(1./(bad*(gamma-1./gamma)))*numgam[i]
            i+=1
            
            
    else:
        print('radtype not supported, try syn or adiab')
        
    
    return numgam

numerical_style = dict(linestyle='--',marker='o',fillstyle='none')

##### CALCULATED ADIABATIC COMPRESSION WITH CONSTANT INJECTION ON A FLAT GAS BACKGROUND #####

# COMPILE CODE
os.chdir('../prob/' + PROBLEM)
call(['python', 'build_adiab_calc.py', '-dir', TMP_DIR])
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

Q0 = 1e5
p = hdr['plaw']
bad = (1./3.)*(-5e-3)
ginjmin = hdr['ginjmin']; ginjmax = hdr['ginjmax']

fig = plt.figure(figsize=(16.18,10))
ax = plt.gca()

# Read in the value of gamma associated with each bin
gammaStart = hdr['NVAR']-hdr['NTEBINS']
gammaString = np.array(hdr['vnams'][gammaStart:])
gammas = gammaString.astype(float)

nDumps = np.size(dfiles)
nSize = np.size(gammas)

nGamma = np.zeros((nDumps,nSize))

for i in range(nDumps):
    dump = io.load_dump(dfiles[i], geom)
    
    # Read in the electron number in each bin in the first zone
    j=0
    for key in gammaString:
        nGamma[i][j] = dump[key][0][0][0]
        j+=1
    
    color=next(ax._get_lines.prop_cycler)['color']
    plt.loglog(gammas,gammas*nGamma[i],label = 't = {:.2e}'.format(float(dump['t'])), **numerical_style, color=color)
    plt.loglog(gammas,gammas*N(float(dump['t']),gammas,'adiab'), color=color)


leg = plt.legend(loc='lower left')
#plt.xlim([4e5,5.5e5])
plt.ylim([1e-11,500])
plt.xlabel('${\gamma}$',fontsize=14)
plt.ylabel('${\gamma}*n({\gamma})$',fontsize=14)
plt.title('Adiabatic Compression on a Flat, Static Background with Scaling (p = 2.5)')

plt.savefig('adiabatic.png', bbox_inches='tight')

# CLEAN UP
# util.safe_remove(TMP_DIR)