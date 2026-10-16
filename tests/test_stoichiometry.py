"""Scientific observation and matched-inventory checks for the POM experiments."""
import unittest
import numpy as np
from splas2 import Config,Physiology,Model
from splas2.model import GPP,RESP_P,UPTAKE_N,INGEST_C,INGEST_N,Z_PROD,RESP_Z
from splas2.stoichiometry import observe,size_metrics
from scripts.run_stoichiometry import registry,build


class StoichiometryObservations(unittest.TestCase):
    def test_observed_processes_reconstruct_independent_rhs(self):
        rng=np.random.default_rng(1729)
        for eta in [.335,.67,1.]:
            m=Model(Config(n_phyto=7,n_zoo=6,depth=1,dilution=.02,inflow_DIN=4),
                    Physiology(assimilation_C=.396,assimilation_N=eta))
            for _ in range(12):
                y=m.initial()*np.exp(rng.normal(0,.8,(m.nvar,1)))
                obs=observe(m,y,10.);f=m.rhs(10.,y)
                pc=m.traits[:,3]*f[m.a,0]+f[m.cr,0]+m.physiology.reserve_N_CN*f[m.nr,0]
                pn=m.traits[:,4]*f[m.a,0]+f[m.nr,0]
                np.testing.assert_allclose(obs['growth_net_C'],pc,rtol=3e-13,atol=1e-14)
                np.testing.assert_allclose(obs['growth_net_N'],pn,rtol=3e-13,atol=1e-14)
                for key,index in [('fixation_C',GPP),('respiration_P_C',RESP_P),('uptake_N',UPTAKE_N),
                                  ('ingestion_C',INGEST_C),('ingestion_N',INGEST_N),('production_Z_C',Z_PROD),('respiration_Z_C',RESP_Z)]:
                    self.assertAlmostEqual(obs[key].sum(),f[m.ntr+index,0],places=11)

    def test_matched_initial_inventories_for_causal_comparisons(self):
        cases=registry();closed=[];minimal=[];community=[]
        for case in cases:
            m=build(case);budget=m.budgets(m.initial(**case['initial']))
            if case['family']=='matched_closed_assembly':closed.append(budget)
            if case['family'].startswith('minimal_'):minimal.append(budget)
            if case['family'].startswith('community_'):community.append(budget)
        for values in [closed,minimal,community]:
            np.testing.assert_allclose(values,np.tile(values[0],(len(values),1)),rtol=1e-14)
        np.testing.assert_allclose(closed[0],[2246.375,8.],rtol=1e-14)

    def test_size_bins_conserve_biomass(self):
        # Equal log-density over the entire fixed physical interval has E=1,
        # including when old bin centres do not coincide with the common grid.
        for n in [35,70,140]:
            x=np.linspace(0,np.log(200),n)
            edges=np.r_[0,.5*(x[:-1]+x[1:]),np.log(200)]
            metrics=size_metrics(np.diff(edges),np.exp(x))
            self.assertAlmostEqual(metrics['evenness'],1.,places=13)

if __name__=='__main__':unittest.main(verbosity=2)
