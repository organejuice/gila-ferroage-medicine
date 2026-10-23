"""Independent toy-network and multiplicity checks for the scientific adapter."""
import importlib.util, unittest
from pathlib import Path
import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('stage030',ROOT/'src/146_network_stage_a_030.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class NetworkAdapterTests(unittest.TestCase):
    def test_multisource_matches_independent_breadth_first_search(self):
        graph=nx.barbell_graph(5,3)
        matrix=nx.to_scipy_sparse_array(graph,nodelist=range(len(graph)),format='csr',dtype=float)
        matrix=csr_matrix((matrix.data,matrix.indices.astype(np.int32),matrix.indptr.astype(np.int32)),shape=matrix.shape)
        seeds=[0,10]
        expected=[min(nx.shortest_path_length(graph,n,s) for s in seeds) for n in graph]
        np.testing.assert_array_equal(m.graph_distance(matrix,seeds),expected)
    def test_sampler_preserves_each_degree_bin_without_duplicates(self):
        graph=nx.barabasi_albert_graph(500,3,seed=27)
        bins,bucket=m.degree_groups(graph,100)
        self.assertEqual(len(bucket),500)
        template=np.array([0,1,2,33,40,77,220,499])
        expected=sorted(bucket[int(n)] for n in template)
        rng=np.random.default_rng(30030)
        for _ in range(100):
            random_set=m.sample_set(template,bins,bucket,rng)
            self.assertEqual(len(set(random_set)),len(template))
            self.assertEqual(sorted(bucket[int(n)] for n in random_set),expected)
    def test_bh_known_hand_calculated_case_and_order(self):
        np.testing.assert_allclose(m.bh([.04,.001,.03,.2]),[.053333333333333,.004,.053333333333333,.2])
        p=np.array([.04,.001,.03,.2]);order=[2,0,3,1]
        np.testing.assert_allclose(m.bh(p[order]),m.bh(p)[order])
    def test_precompiled_draws_are_bitwise_identical(self):
        graph=nx.barabasi_albert_graph(500,3,seed=27)
        bins,bucket=m.degree_groups(graph,100);template=np.array([0,1,2,33,40,77,220,499])
        original=np.random.default_rng(30030);optimized=np.random.default_rng(30030)
        plan=m.sample_plan(template,bucket)
        for _ in range(100):
            np.testing.assert_array_equal(m.sample_set(template,bins,bucket,original),m.draw_plan(plan,bins,optimized))
    def test_unreachable_nodes_and_self_distance(self):
        graph=nx.Graph([(0,1),(2,3)])
        matrix=nx.to_scipy_sparse_array(graph,nodelist=range(4),dtype=float,format='csr')
        matrix=csr_matrix((matrix.data,matrix.indices.astype(np.int32),matrix.indptr.astype(np.int32)),shape=matrix.shape)
        d=m.graph_distance(matrix,[0]);self.assertEqual(d[0],0);self.assertEqual(d[1],1);self.assertTrue(np.isinf(d[2]))

if __name__=='__main__':unittest.main()
