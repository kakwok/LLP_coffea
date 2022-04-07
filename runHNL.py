from coffea import hist, processor
from optparse import OptionParser
import coffea
import pickle,glob
import time
from coffea.nanoevents import NanoEventsFactory, BaseSchema

def runLocal(outf="test.pickle",fileset="test.json",isElectronChannel=True,full=False,saveSkim=False,debug=False):
    #from HNLprocessor.HNLproc_3 import MyProcessor
    from HNLprocessor.HNLproc_4 import MyProcessor

    if full:      
        out = processor.run_uproot_job(
            fileset,
            treename="MuonSystem",
            processor_instance=MyProcessor(isElectronChannel,saveSkim,debug),
            executor=processor.iterative_executor,
            executor_args={
                "schema": BaseSchema,
            },
            chunksize=10000,
        )
    else:
        out = processor.run_uproot_job(
            fileset,
            treename="MuonSystem",
            processor_instance=MyProcessor(isElectronChannel,saveSkim,debug),
            executor=processor.iterative_executor,
            executor_args={
                "schema": BaseSchema,
            },
            maxchunks = 2,
            chunksize=100,
        )

    print(out)
    with open(outf,'wb') as f:
        pickle.dump(out,f)
    return

def runLPC(outf="test.pickle",fileset="test.json",isElectronChannel=True,nJobs=4,saveSkim=False):
    import time
    from distributed import Client
    from lpcjobqueue import LPCCondorCluster
    tic = time.time()
    #cluster = LPCCondorCluster(log_directory="/uscms/home/kkwok/log")
    cluster = LPCCondorCluster()

    # minimum > 0: https://github.com/CoffeaTeam/coffea/issues/465
    cluster.adapt(minimum=4, maximum=100)
    client = Client(cluster)

    exe_args = {
        "client": client,
        "savemetrics": True,
        "schema": BaseSchema,
        "align_clusters": True,
    }

    client.upload_file("HNLprocessor.zip")


    #from HNLprocessor.HNLproc_3 import MyProcessor
    from HNLprocessor.HNLproc_4 import MyProcessor

    print("Waiting for at least one worker...")
    client.wait_for_workers(4)
    hists, metrics = processor.run_uproot_job(
        fileset,
        treename="MuonSystem",
        processor_instance=MyProcessor(isElectronChannel,saveSkim),
        executor=processor.dask_executor,
        executor_args=exe_args,
        chunksize=10000,
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

    parser = OptionParser()
    parser.add_option('--test', dest='test', action='store_true',default = False, help='Run local test with small fileset')
    parser.add_option('--local', dest='local', action='store_true',default = False, help='Run local test with 1 chunk of full fileset')
    parser.add_option('--condor', dest='condor', action='store_true',default = False, help='Run local test with 1 chunk of full fileset')
    parser.add_option('--muon', dest='muon', action='store_true',default = False, help='Run muon channel')
    parser.add_option('--full', dest='full', action='store_true',default = False, help='Run full file chunks')
    parser.add_option('--saveSkim', dest='saveSkim', action='store_true',default = False, help='Save skim selections')
    parser.add_option('--debug', dest='debug', action='store_true',default = False, help='run with debug')
    parser.add_option('--fileset', dest='fileset', default = "test.json", help='input file json')
    parser.add_option('--nJobs', dest='nJobs', default = 4, type=int, help='number of workers in condor')
    parser.add_option('-o', dest='outf', default='HNL_histograms.pickle', help='collection of histograms')

    (options, args) = parser.parse_args()
    outf    = options.outf 
    saveSkim    = options.saveSkim
    fileset = options.fileset
    isElectronChannel = not options.muon

    print(" Fileset = ", fileset)
    print(" isElectronChannel = ", isElectronChannel)
    print(" outf              = ", outf)
    print(" saveSkim          = ", saveSkim)

    if options.test:
        runLocal("test.pickle","test.json",isElectronChannel)
    elif options.local:
        print("full               = ", options.full)
        runLocal(outf,fileset,isElectronChannel,options.full,saveSkim,options.debug)
    elif options.condor:
        print(" Using nJobs       = ", options.nJobs)
        runLPC(outf,fileset,isElectronChannel,options.nJobs,saveSkim)
