from objects.complex import *

complex_ids = os.listdir(BASE_DATA_PATH)
reprocess = False
n_of_patchdock_results = 200

for complex_id in complex_ids:
    print("Processing %s" % complex_id)
    BenchmarkComplex(complex_id, ComplexType.zdock_benchmark_unbound, reprocess)
    BenchmarkComplex(complex_id, ComplexType.zdock_benchmark_bound, reprocess)
    """
    for i in range(n_of_patchdock_results):
        PatchDockComplex(complex_id, i + 1, reprocess)
    """