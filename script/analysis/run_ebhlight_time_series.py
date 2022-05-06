from __future__ import print_function
import numpy as np
import h5py
import glob
import sys

# ebhlight 1D time series of L_cool = \int dV \Lambda, Lbol integrated over \Omega, total # of photons emitted and on grid

sim_name = sys.argv[1]
print('ebhlight time series: ',sim_name)
fname='ebhlight_timeseries_'+sim_name
files=sorted(glob.glob('dumps/dump*.h5'))

gf = h5py.File('dumps/grid.h5','r')
gdet = gf['gdet'][()]
r = gf['Xbl'][()][:,0,0,1]
gf.close()

t = []; lcool = []; lbol = []; nsph = []; nem = []
for f in files:
    print(f)
    f = h5py.File(f,'r')
    Lam=-f['Jrad'][()][0]*f['dt'][()]/f['dump_cadence'][()]
    Rmunu=f['Rmunu'][()]
    dx1=f['geom/dx1'][()]; dx2=f['geom/dx2'][()]; dx3=f['geom/dx3'][()]
    M_unit=f['/header/units/M_unit'][()]
    L_unit=f['header/units/L_unit'][()]
    T_unit=f['/header/units/T_unit'][()]
    t.append(f['t'][()])
    lumunit=M_unit*L_unit**2./T_unit**3.
    Nsph=f['Nsph'][()]; Nem=f['Nem'][()]
    lcool.append(np.sum(gdet*Lam)*dx1*dx2*dx3*lumunit)
    Lr = np.sum(np.sum(-gdet*f['Rmunu'][:,:,:,1,0],-1),-1)*dx2*dx3*lumunit
    rad_indx=np.max(np.where(r < f['geom/mmks/r_out_rad'][()]))
    lbol.append(np.mean(Lr[rad_indx-10:rad_indx]))
    nsph.append(np.sum(Nsph)); nem.append(np.sum(Nem))
    f.close()

np.save(fname,(t,lcool,lbol,nsph,nem))
