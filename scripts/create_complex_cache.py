from Constants import TRAIN_COMPLEX_IDS
from objects.complex import *

complex_ids = TRAIN_COMPLEX_IDS
recache = False
n_of_patchdock_results = 200

for complex_id in complex_ids:
    print("Caching %s" % complex_id)
    BenchmarkComplex(complex_id, ComplexType.zdock_benchmark_unbound, recache)
    BenchmarkComplex(complex_id, ComplexType.zdock_benchmark_bound, recache)
    for i in range(n_of_patchdock_results):
        PatchDockComplex(complex_id, i + 1, recache)
