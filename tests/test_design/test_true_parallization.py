from dasi.design import Design, LibraryDesign
from multiprocessing import Pool
import pytest

@pytest.mark.parametrize('i', [1, 2, None])
def test_compile_with_filtered_keys(i):
    design = Design.fake(3)
    design._blast()
    query_keys = design.query_keys
    design.assemble_graphs(n_jobs=1, query_keys=query_keys[:i])
    design.pre_process_graphs(**{})
    design.optimize()
    assert len(design.results) == i

def run(args):
    design, qk = args
    design.assemble_graphs(n_jobs=1, query_keys=qk)
    design.pre_process_graphs(**{})
    design.optimize(n_jobs=1)
    return design.results

def test_true_parallelization():
    design = LibraryDesign.fake(10)
    design.DEFAULT_N_JOBS = 1
    design._blast()
    design._share_query_blast()
    design._expand_from_synthesized()
    design.update_library_metadata()
    query_keys = design.query_keys
    chunks = []
    for i in range(10):
        chunks.append(query_keys[i:i+1])
    # chunks = [query_keys]
    args = list(zip([design] * len(chunks), chunks))
    with Pool(processes=min(len(chunks), 16)) as pool:
        results = pool.map(run, args)
    print(results)