#from HNLproc import MyProcessor
#from HNLproc_2 import MyProcessor
#from HNLproc_3 import MyProcessor
from coffea import hist, processor
import coffea
import pickle,glob
import time

def runLocal(outf="test.pickle",fileset="test.json"):
    from HNLprocessor.HNLproc_3 import MyProcessor
    out = processor.run_uproot_job(
        #"samplefiles.json",
        #"all_samples.json",
        #"signals_skim.json",
        fileset,
        treename="MuonSystem",
        processor_instance=MyProcessor(),
        executor=processor.iterative_executor,
        executor_args={
            "schema": None,
        },
        maxchunks=2,
        chunksize=100
    )
    with open(outf,'wb') as f:
        pickle.dump(out,f)
    return

def runLPC(outf="test.pickle",fileset="test.json"):
    import time
    from distributed import Client
    from lpcjobqueue import LPCCondorCluster
    tic = time.time()
    cluster = LPCCondorCluster()
    # minimum > 0: https://github.com/CoffeaTeam/coffea/issues/465
    cluster.adapt(minimum=4, maximum=10)
    client = Client(cluster)

    exe_args = {
        "client": client,
        "savemetrics": True,
        "schema": None,
        "align_clusters": True,
    }

    client.upload_file("HNLprocessor.zip")


    from HNLprocessor.HNLproc_3 import MyProcessor

    print("Waiting for at least one worker...")
    client.wait_for_workers(4)
    hists, metrics = processor.run_uproot_job(
        fileset,
        treename="MuonSystem",
        processor_instance=MyProcessor(),
        executor=processor.dask_executor,
        executor_args=exe_args,
        # remove this to run on the whole fileset:
        #maxchunks=10,
    )


    elapsed = time.time() - tic
    print(f"Output: {hists}")
    print(f"Metrics: {metrics}")
    print(f"Finished in {elapsed:.1f}s")
    print(f"Events/s: {metrics['entries'] / elapsed:.0f}")

    with open(outf,'wb') as f:
        pickle.dump(hists,f)
    return

if __name__ == '__main__':

    #outf = "HNL_histograms_signals_Jul28.pickle"
    #outf = "HNL_histograms_signals_Aug5.pickle"
    outf = "HNL_histograms_all_Sep17.pickle"
    fileset = "signals_skim.json"
    #fileset = "test.json"

    runLPC(outf,fileset)
    #runLocal(outf,fileset)
    #runLocal("test.pickle","test.json")
