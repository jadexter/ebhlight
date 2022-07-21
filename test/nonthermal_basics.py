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

##### SYNCHROTRON COOLING WITH CONSTANT INJECTION ON A FLAT GAS BACKGROUND #####

# COMPILE CODE
os.chdir('../prob/' + PROBLEM)
call(['python', 'build_synchrotron.py', '-dir', TMP_DIR])
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
bs = 1.292e-11*(200**2)
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
        nGamma[i][j] = dump[key][0]
        j+=1
    
    
    color=next(ax._get_lines.prop_cycler)['color']    
    plt.loglog(gammas,gammas*nGamma[i],label = 't = {:.2e}'.format(float(dump['t'])),color=color, **numerical_style)
    plt.loglog(gammas,gammas*N(float(dump['t']),gammas,'syn'),color=color)
    # plt.axvline(pow(((1e-5)+((1.292e-11)*pow(200,2)*float(dump['t']))),-1))

plt.loglog(gammas[25:35], (3e11*gammas**(-3.5))[25:35],  color='r',ls='--')
plt.loglog(gammas[30:40], (3e3*gammas**(-2.5))[-22:-12], color='r',ls='--')
plt.loglog(gammas[9:15],  (3e4*gammas**(-1))[9:15],      color='r',ls='--')

plt.text(gammas[30], 0.5,   '${\gamma}^{-3.5}$')
plt.text(gammas[35], 5e-10, '${\gamma}^{-2.5}$')
plt.text(gammas[12], 300,   '${\gamma}^{-1}$')

leg = plt.legend(loc='upper right')
#plt.xlim([4e5,5.5e5])
plt.ylim([1e-12,1.5e3])
plt.xlabel('${\gamma}$',fontsize=14)
plt.ylabel('${\gamma}*n({\gamma})$',fontsize=14)
plt.title('Synchrotron Cooling on a Flat, Static Background (p = 3.5)')

plt.savefig('synchrotron.png', bbox_inches='tight')

# Only run this mkdir piece the first time
util.safe_remove('dump_storage')
call(['mkdir', './dump_storage'])

call(['mv', './' + TMP_DIR + '/dumps','./dump_storage'])
call(['mv', './dump_storage/dumps', './dump_storage/synchrotron_dumps'])
# CLEAN UP
util.safe_remove(TMP_DIR)



##### ADIABATIC COMPRESSION WITH CONSTANT INJECTION ON A FLAT GAS BACKGROUND #####

# COMPILE CODE
os.chdir('../prob/' + PROBLEM)
call(['python', 'build_adiab.py', '-dir', TMP_DIR])
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
        nGamma[i][j] = dump[key][0]
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
call(['mv', './' + TMP_DIR + '/dumps','./dump_storage'])
call(['mv', './dump_storage/dumps', './dump_storage/adiabatic_dumps'])
util.safe_remove(TMP_DIR)



##### SIMPLE TESTING OF SPATIAL EVOLUTION IN A SHOCKTUBE #####

# COMPILE CODE
os.chdir('../prob/' + PROBLEM)
call(['python', 'build_flux.py', '-dir', TMP_DIR])
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

gammaStart = hdr['NVAR']-hdr['NTEBINS']
gammaString = np.array(hdr['vnams'][gammaStart:])
gammas = gammaString.astype(float)

tscale = 1.e-2
gam = 1.4

rho_code = dump['RHO'][:,0,0]
P_code   = (gam-1.)*dump['UU'][:,0,0]/(tscale*tscale)
u_code   = dump['U1'][:,0,0]/(tscale)
#N1 = len(rho)
x_code = geom['x'][:,0,0]
n_code = dump[gammaString[16]][:,0,0]

# GET ANALYTIC SOLUTION (ADAPTED FROM BRUCE FRYXELL'S exact_riemann.f)
x0 = 0.5
t = 0.25
rho1 = 1.; P1 = 1.; u1 = 0.
rho5 = 0.125; P5 = 0.1; u5 = 0.
cs1 = np.sqrt(gam*P1/rho1)
cs5 = np.sqrt(gam*P5/rho5)
Gam = (gam-1.)/(gam+1.)
beta = (gam-1.)/(2.*gam)
def func(x):
  P3 = x[0]
  u4 = (P3-P5)*np.sqrt((1.-Gam)/(rho5*(P3 + Gam*P5)))
  u2 = (P1**beta - P3**beta)*np.sqrt((1.-Gam**2.)*P1**(1./gam)/(Gam**2*rho1))
  return u2-u4
P3 = root(func, [(P1+P5)/2.]).x[0]
P4 = P3
rho3 = rho1*(P3/P1)**(1./gam)
rho4 = rho5*(P4 + Gam*P5)/(P5 + Gam*P4)
u4 = cs5*(P4/P5-1)/(gam*np.sqrt(1. + (gam+1.)/(2.*gam)*(P4/P5-1.)))
ushock = cs5*np.sqrt(1. + (gam+1.)/(2.*gam)*(P4/P5-1.))
u3 = u4
cs3 = np.sqrt(gam*P3/rho3)
xsh = x0 + ushock*t
xcd = x0 + u3*t
xft = 0.5 + (u3-cs3)*t
xhd = 0.5 - cs1*t
N = 1024
x = np.linspace(0, 1, 1024)
rho = np.zeros(N)
P   = np.zeros(N)
u   = np.zeros(N)
for n in range(N):
  if x[n] < xhd:
    rho[n] = rho1
    P[n]   = P1
    u[n]   = u1
  elif x[n] < xft:
    u[n]   = 2./(gam+1.)*(cs1 + (x[n] - x0)/t)
    fac    = 1. - 0.5*(gam-1.)*u[n]/cs1
    rho[n] = rho1*fac**(2./(gam-1.))
    P[n]   = P1*fac**(2.*gam/(gam-1.))
  elif x[n] < xcd:
    rho[n] = rho3
    P[n]   = P3
    u[n]   = u3
  elif x[n] < xsh:
    rho[n] = rho4
    P[n]   = P4
    u[n]   = u4
  else:
    rho[n] = rho5
    P[n]   = P5
    u[n]   = u5

code_col = 'r'; code_ls = ''; code_mrk = '.'
sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(1,2,1)
ax.plot(x_code, rho_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(x, rho, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Density')
plt.ylim([0, 1.1])

ax = fig.add_subplot(1,2,2)
ax.plot(x_code, n_code, color=code_col, linestyle=code_ls, marker=code_mrk)
plt.xlabel('x'); plt.ylabel('N(gamma)')

plt.subplots_adjust(wspace=0.15)

plt.savefig('fluxtest.png', bbox_inches='tight')


# CLEAN UP
call(['mv', './' + TMP_DIR + '/dumps','./dump_storage'])
call(['mv', './dump_storage/dumps', './dump_storage/flux_dumps'])
util.safe_remove(TMP_DIR)

util.safe_remove('test_images')
call(['mkdir', './test_images'])
call(['mv', './synchrotron.png', './test_images'])
call(['mv', './adiabatic.png',   './test_images'])
call(['mv', './fluxtest.png',    './test_images'])
