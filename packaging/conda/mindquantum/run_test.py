import numpy as np
from scipy import sparse
from mindquantum.algorithm.qaia import CFC

solver = CFC(sparse.csr_matrix([[0., 1.], [1., 0.]]),
             x=np.full((2, 4), 0.001), n_iter=20, batch_size=4, backend='cpu-float32')
solver.update()
assert solver.backend == 'cpu-float32'
assert solver.x.shape == (2, 4)
assert np.isfinite(solver.x).all()
print('MindQuantum CPU CFC passed')
