"""Checks on scientific observables, not acceptance tests for ecological patterns."""
import unittest
import numpy as np
from splas2 import Config,Physiology,Model
from splas2.pom_observation import nitrogen_network,rebin
from scripts.run_pom_validation import registry,build


class POMObservation(unittest.TestCase):
    def test_nitrogen_routing_matches_independent_rhs(self):
        rng=np.random.default_rng(3519)
        for closure in [0.,1.]:
            for eta in [.335,.67,1.]:
                p=Physiology(zoo_closure=closure,assimilation_N=eta,mortality_recycled=.6)
                m=Model(Config(n_phyto=7,n_zoo=6,depth=1,dilution=.02,inflow_DIN=4),p)
                for _ in range(5):
                    y=m.initial()*np.exp(rng.normal(0,.6,(m.nvar,1)))
                    f=m.rhs(0.,y);o=nitrogen_network(m,y)
                    expected=np.array([f[0,0],(m.traits[:,4]*f[m.a,0]+f[m.nr,0]).sum(),
                        f[m.z,0].sum()/p.zoo_CN,f[m.dn,0].sum(),f[3,0]])
                    np.testing.assert_allclose(o['derivative'],expected,rtol=8e-13,atol=1e-13)

    def test_new_ensemble_has_matched_initial_inventories(self):
        budgets=[]
        for c in registry():
            m=build(c);budgets.append(m.budgets(m.initial(**c['initial'])))
        np.testing.assert_allclose(budgets,np.tile([2205.0225,1.77],(88,1)),rtol=2e-14)

    def test_rebin_conserves_both_trophic_size_domains(self):
        for lo,hi in [(1.,200.),(2.,1000.)]:
            for n in [35,70,140]:
                sizes=np.geomspace(lo,hi,n);b=np.arange(1,n+1,dtype=float)
                for nb in [7,14,28]:
                    _,q=rebin(b,sizes,lo,hi,nb)
                    self.assertAlmostEqual(q.sum()/b.sum(),1.,places=14)

if __name__=='__main__':unittest.main(verbosity=2)
