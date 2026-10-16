"""Historical didactic Monod NPZ, with a conservative fixed C:N/DIC ledger.

Not a reproduction of the old Fortran executable. No quota acclimation, no
maintenance, no detrital delay: an explicit competing simplified hypothesis.
"""
import numpy as np
from .model import njit

@njit(cache=True)
def fixed_rhs(y,x,xz,phi,mort=.1,sink=.0001,eps=.33,xi=.2,ratio=6.625):
 n=len(x);nz=len(xz);dy=np.zeros_like(y);loss=np.zeros(n)
 for j in range(nz):
  av=0.
  for i in range(n):av+=phi[i,j]*y[2+i]
  f=26*xz[j]**(-.4)*y[2+n+j]/(3.+av)
  ing=0.
  for i in range(n):loss[i]+=f*phi[i,j];ing+=f*phi[i,j]*y[2+i]
  zm=xi*y[2+n+j]*y[2+n:].sum()
  dy[2+n+j]=eps*ing-zm
  recycled=(1-eps)*ing+zm;dy[0]+=recycled;dy[1]+=ratio*recycled
 for i in range(n):
  grow=2.6*x[i]**(-.45)*y[0]/(.1*x[i]+y[0])*y[1]/(10+y[1])*y[2+i]
  death=(mort+sink*x[i]**1.17)*y[2+i]
  dy[2+i]=grow-loss[i]*y[2+i]-death
  dy[0]+=-grow+death;dy[1]+=ratio*(-grow+death)
 return dy

@njit(cache=True)
def _fixed_advance(y,x,xz,phi,end,maxstep,mort,sink):
 t=0.;steps=0
 while t<end-1e-12:
  f=fixed_rhs(y,x,xz,phi,mort,sink);h=min(maxstep,end-t)
  # Gross donor-turnover bound, not just net dDIN: prevents a positive
  # oscillatory numerical equilibrium in rapidly recycled nutrients.
  gross=0.
  for i in range(len(x)):
   gross+=2.6*x[i]**(-.45)*y[0]/(.1*x[i]+y[0])*y[1]/(10+y[1])*y[2+i]
  if gross>0:h=min(h,.8*y[0]/gross,.8*y[1]/(6.625*gross))
  for i in range(len(y)):
   if f[i]<0:h=min(h,.85*y[i]/(-f[i]))
  if h<1e-12:raise RuntimeError('Fixed-NPZ positivity underflow')
  for _ in range(40):
   a=y+h*f;f2=fixed_rhs(a,x,xz,phi,mort,sink);b=a+h*f2
   if np.all(b>=0):break
   h*=.5
  y=.5*(y+b);t+=h;steps+=1
 return y,steps

def run_fixed(n_zoo=0,days=3652.,max_step=.04,mortality=.1,sinking=.0001):
 if not 0<=n_zoo<=70 or days<=0 or max_step<=0:raise ValueError('Invalid fixed-NPZ experiment')
 x=np.geomspace(1,200,70);xz=np.geomspace(2,1000,70)[:n_zoo]
 opt=.16*xz*np.exp(-.02*np.log(xz)**2)
 phi=np.exp(-((np.log10(x[:,None])-np.log10(opt[None,:]))/.25)**2)
 y=np.r_[1.,4000.,np.full(70,.1),np.full(n_zoo,.01)];initial=y.copy()
 y,steps=_fixed_advance(y,x,xz,phi,days,max_step,mortality,sinking)
 return {'initial':initial,'final':y,'sizes':x,'zsizes':xz,'steps':steps,
   'mass_error_N':float(y[0]+y[2:].sum()-initial[0]-initial[2:].sum()),
   'mass_error_C':float(y[1]+6.625*y[2:].sum()-initial[1]-6.625*initial[2:].sum())}
