"""Independent element accounting, analytical transport and falsifying edge cases."""
import unittest
from dataclasses import replace
import numpy as np
from splas2 import Model,Config,Physiology,simulate
from splas2.model import *

class ScientificInvariants(unittest.TestCase):
 def test_reference_allometry_and_units(self):
  self.assertAlmostEqual(float(montagnes_carbon_pg(100)),.109*100**.991)
  m=Model(Config(n_phyto=1,n_zoo=0,phyto_min=10,phyto_max=10))
  self.assertAlmostEqual(m.traits[0,1]*120,float(montagnes_carbon_pg(np.pi/6*1000)))
  self.assertAlmostEqual(float(m.Dface[0]),0.)
  self.assertAlmostEqual(Model(Config(diffusivity_m2_s=1e-4)).Dface[0],8.64)
 def test_inventory_independent_of_resolution(self):
  inventories=[]
  for n in [7,35,70,140]:
   m=Model(Config(n_phyto=n,n_zoo=n));inventories.append(m.budgets(m.initial()))
  np.testing.assert_allclose(inventories,np.repeat([inventories[0]],4,axis=0),rtol=1e-14)
 def test_random_local_and_spatial_element_balance(self):
  rng=np.random.default_rng(1772)
  for open_case in [False,True]:
   m=Model(Config(n_phyto=9,n_zoo=7,n_layers=5,depth=25,diffusivity_m2_s=1e-4,deep_diffusivity_m2_s=1e-3,
      export_bottom=open_case,nutrient_supply=.2 if open_case else 0,bottom_N_flux=4 if open_case else 0,
      dilution=.03 if open_case else 0,dic_exchange=.1 if open_case else 0),Physiology(mortality_recycled=.3 if open_case else 1.))
   for _ in range(100):
    y=m.initial()*np.exp(rng.normal(0,1,(m.nvar,5)))
    y[m.dc]=rng.uniform(0,2,(m.nD,5));y[m.dn]=rng.uniform(0,1,(m.nD,5));y[DOC:DON+1]=rng.uniform(0,2,(2,5))
    dy=m.rhs(float(rng.uniform(0,400)),y)
    # Element sums from independently reconstructed biological inventories.
    ca=(m.traits[:,3,None]*dy[m.a]+dy[m.cr]+m.physiology.reserve_N_CN*dy[m.nr]).sum()+dy[DIC].sum()+dy[DOC].sum()+dy[m.z].sum()+dy[m.dc].sum()
    na=(m.traits[:,4,None]*dy[m.a]+dy[m.nr]).sum()+dy[DIN].sum()+dy[DON].sum()+dy[m.z].sum()/m.physiology.zoo_CN+dy[m.dn].sum()
    d=dy[m.ntr:].sum(axis=1)
    self.assertLess(abs(ca+d[EXPORT_C]+d[HTL_C]+d[OUT_C]-d[INPUT_C]),1e-9)
    self.assertLess(abs(na+d[EXPORT_N]+d[HTL_N]+d[OUT_N]-d[INPUT_N]),1e-10)
 def test_positive_orthant_boundary(self):
  m=Model(Config(n_phyto=5,n_zoo=4,n_layers=3,depth=15,export_bottom=True,diffusivity_m2_s=1e-4))
  for r in range(m.ntr):
   y=m.initial();y[r,1]=0
   if m.a.start<=r<m.a.stop:y[m.cr.start+r-m.a.start,1]=0;y[m.nr.start+r-m.a.start,1]=0
   self.assertGreaterEqual(m.rhs(0,y)[r,1],-1e-13)
 def test_no_energetically_free_nitrogen_uptake(self):
  m=Model(Config(n_phyto=1,n_zoo=0,surface_PAR=0))
  y=m.initial();y[m.cr]=0
  f=m.rhs(0,y)
  self.assertEqual(f[m.ntr+UPTAKE_N,0],0)
  self.assertGreater(f[m.ntr+RESP_P,0],0)
  self.assertLess(f[m.a,0].item(),0)
  r=simulate(m,10,y0=y,save_interval=1)
  self.assertLess(r.diagnostics()['max_relative_mass_error_C'],1e-12)
  self.assertGreaterEqual(r.state[:,:m.ntr].min(),0)
 def test_nitrogen_free_reserves_can_support_initial_division(self):
  m=Model(Config(n_phyto=1,n_zoo=0));y=m.initial(DIN_value=0)
  f=m.rhs(0,y);self.assertEqual(f[m.ntr+UPTAKE_N,0],0)
  self.assertGreater(f[m.a,0].item()+m.physiology.mortality*y[m.a,0].item(),0)
 def test_no_new_organic_c_without_dic_or_reserves(self):
  m=Model(Config(n_phyto=1,n_zoo=0));y=m.initial(DIC_value=0);y[m.cr]=0
  f=m.rhs(0,y);self.assertEqual(f[m.ntr+GPP,0],0)
  self.assertLessEqual(f[m.a,0].item(),0)
 def test_ingestion_bounded_and_nonsymmetric_kernel(self):
  m=Model(Config(n_phyto=8,n_zoo=8));y=m.initial(phyto_N=1e6,zoo_N=.1)
  f=m.rhs(0,y)
  maxing=(26*m.zsizes**(-.4)*y[m.z,0]).sum()
  self.assertLessEqual(f[m.ntr+INGEST_C,0],maxing*(1+1e-12))
  bc,_=m.biomass(y);available=m.phi.T@bc[:,0]
  oracle=(26*m.zsizes**(-.4)*y[m.z,0]*available/(m.config.grazing_half_C+available)).sum()
  self.assertAlmostEqual(f[m.ntr+INGEST_C,0],oracle,places=12)
 def test_homogeneous_0d_light_and_exact_annual_dose(self):
  m=Model(Config(n_phyto=1,n_zoo=0,surface_PAR=85,forcing=4));y=m.initial()
  self.assertEqual(m.light(0,y)[0][0],85);self.assertEqual(m.light(.75,y)[0][0],0)
  t=np.linspace(0,365,7301)
  for force in [1,2,3]:
   v=np.array([surface_light(z,100,force) for z in t]);self.assertAlmostEqual(np.trapezoid(v,t)/365,75,places=8)
 def test_diffusion_analytic_mode_and_variable_D_conservation(self):
  p=Physiology(mu_scale=0,photo_scale=0,uptake_scale=0,mortality=0,maintenance=0,exud_C=0,exud_N=0)
  m=Model(Config(n_phyto=1,n_zoo=0,n_layers=20,depth=20,diffusivity_m2_s=1e-4),p)
  y=m.initial(phyto_N=0,zoo_N=0);mode=np.cos(np.pi*m.depth/20);y[DIN]=2+mode
  r=simulate(m,2,y0=y,max_step=.01,times=[0,2])
  eigen=4*8.64/m.dz**2*np.sin(np.pi/40)**2
  np.testing.assert_allclose(r.state[-1,DIN],2+mode*np.exp(-eigen*2),rtol=3e-7,atol=3e-7)
  self.assertLess(r.diagnostics()['max_relative_mass_error_N'],1e-12)
 def test_bottom_export_matches_exact_initial_flux(self):
  m=Model(Config(n_phyto=3,n_zoo=0,n_layers=4,depth=20,export_bottom=True))
  y=m.initial();f=m.rhs(0,y)
  bc,bn=m.biomass(y)
  expected=(m.w[m.a]*bc[:,-1]).sum()/m.dz
  self.assertAlmostEqual(f[m.ntr+EXPORT_C,-1],expected,places=12)
 def test_open_integrated_budgets_and_time_convergence(self):
  m=Model(Config(n_phyto=7,n_zoo=5,n_layers=4,depth=20,export_bottom=True,diffusivity_m2_s=1e-4,
       nutrient_supply=.1,bottom_N_flux=3,dilution=.02,dic_exchange=.1),Physiology(mortality_recycled=.3))
  from scipy.integrate import solve_ivp
  y0=m.initial();shape=y0.shape
  oracle=solve_ivp(lambda t,y:m.rhs(t,y.reshape(shape)).ravel(),[0,3],y0.ravel(),method='DOP853',rtol=1e-11,atol=1e-13)
  self.assertTrue(oracle.success)
  target=oracle.y[:,-1].reshape(shape)[:m.ntr]
  rr=[simulate(m,3,max_step=.05,times=[0,3],rtol=tol,atol=tol*1e-4) for tol in [1e-3,1e-4,1e-5]]
  for r in rr:
   self.assertLess(r.diagnostics()['max_relative_mass_error_C'],1e-11)
   self.assertLess(r.diagnostics()['max_relative_mass_error_N'],1e-11)
   self.assertGreaterEqual(r.state[:,:m.ntr].min(),0)
  errors=[np.linalg.norm(r.state[-1,:m.ntr]-target) for r in rr]
  self.assertLess(errors[-1],errors[0]/3)
  self.assertLess(errors[-1]/np.linalg.norm(target),1e-6)
 def test_Rstar_equilibrium_against_independent_algebra(self):
  from splas2.npz_fixed import run_fixed
  from scipy.optimize import brentq
  r=run_fixed(0);nt=8.;ct=4046.375
  root=brentq(lambda N:2.6*N/(.1+N)*(ct-6.625*(nt-N))/(ct-6.625*(nt-N)+10)-.1001,1e-8,1)
  self.assertLess(abs(r['final'][0]-root),1e-8)
  self.assertLess(abs(r['mass_error_N']),1e-8)
 def test_organic_N_reserve_carries_C(self):
  m=Model(Config(n_phyto=6,n_zoo=0));rng=np.random.default_rng(94)
  for _ in range(100):
   y=m.initial()*np.exp(rng.normal(0,2,(m.nvar,1)));bc,bn=m.biomass(y)
   self.assertGreaterEqual((bc-m.physiology.reserve_N_CN*bn).min(),-1e-10)
  m=Model(Config(n_phyto=1,n_zoo=0,phyto_min=5,phyto_max=5,surface_PAR=0))
  r=simulate(m,10,y0=m.initial(DIN_value=100),times=[0,10])
  bc,bn=m.biomass(r.state[-1]);self.assertGreaterEqual((bc/bn).min(),4)
 def test_low_light_slope_not_double_N_penalized(self):
  m=Model(Config(n_phyto=1,n_zoo=0,phyto_min=5,phyto_max=5,surface_PAR=1e-6))
  y=m.initial();y[m.nr]*=.1
  a=y[m.a,0].item();cr=y[m.cr,0].item();nr=y[m.nr,0].item();p=m.physiology;t=m.traits[0]
  status=(t[4]*a+nr)/(t[2]*a);theta=p.theta_min+(p.theta_max-p.theta_min)*status/(1+1e-6/p.acclimation_PAR)
  capacity=max(0,1-(cr/(t[8]*a-t[3]*a-p.reserve_N_CN*nr))**2)
  expected=a*t[1]*p.alpha_Chl*theta*2200/(2200+p.half_DIC)*capacity
  self.assertAlmostEqual(m.rhs(0,y)[m.ntr+GPP,0]/1e-6,expected,places=6)
 def test_discontinuous_day_night_refinement(self):
  m=Model(Config(n_phyto=1,n_zoo=0,forcing=4,surface_PAR=85))
  a=simulate(m,2.3,max_step=.037,times=[0,2.3]);b=simulate(m,2.3,max_step=.004,times=[0,2.3])
  np.testing.assert_allclose(a.state[-1,:m.ntr],b.state[-1,:m.ntr],rtol=.002,atol=3e-5)
 def test_implicit_uptake_converges_to_full_explicit_reference(self):
  m=Model(Config(n_phyto=1,n_zoo=0,phyto_min=5,phyto_max=5))
  ref=simulate(m,5,max_step=.001,times=[0,5],rtol=1e-7,atol=1e-11)
  rr=[simulate(m,5,max_step=h,times=[0,5],method='imex') for h in [.05,.025,.0125]]
  errors=[np.linalg.norm(r.state[-1,:m.ntr]-ref.state[-1,:m.ntr]) for r in rr]
  self.assertLess(errors[-1],errors[0]/2.5)
  for r in rr:
   self.assertGreaterEqual(r.state[:,:m.ntr].min(),0.)
   self.assertLess(r.diagnostics()['max_relative_mass_error_N'],1e-11)
 def test_implicit_energetic_cap_and_extinction(self):
  for abundance_scale in [1.,1e-200]:
   m=Model(Config(n_phyto=4,n_zoo=3,n_layers=3,depth=15,nutrient_supply=1,diffusivity_m2_s=1e-4))
   y=m.initial(DIN_value=100);y[m.a]*=abundance_scale;y[m.cr]*=abundance_scale;y[m.nr]*=abundance_scale
   r=simulate(m,2,y0=y,method='imex',save_interval=.2)
   self.assertGreaterEqual(r.state[:,:m.ntr].min(),0.)
   self.assertLess(r.diagnostics()['max_relative_mass_error_C'],1e-11)
   self.assertLess(r.diagnostics()['max_relative_mass_error_N'],1e-11)
 def test_particle_rebinning_preserves_both_elements(self):
  for n in [1,2,4,8]:
   m=Model(Config(n_detritus=n));np.testing.assert_allclose(m.pw.sum(axis=1),1,atol=1e-15);np.testing.assert_allclose(m.zw.sum(axis=1),1,atol=1e-15)
 def test_validate_failures(self):
  for cfg in [Config(depth=0),Config(diffusivity_m2_s=-1),Config(n_phyto=0)]:
   with self.assertRaises(ValueError):Model(cfg)
  with self.assertRaises(ValueError):Model(physiology=Physiology(acclimation_PAR=0))
if __name__=='__main__':unittest.main(verbosity=2)
