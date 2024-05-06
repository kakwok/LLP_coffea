from coffea import hist, processor
from optparse import OptionParser
import coffea
import pickle,glob
import time
from coffea.nanoevents import NanoEventsFactory, BaseSchema
from HNLprocessor.HNLproc_4 import MyProcessor
from HNLprocessor.TTbar_proc import ttbarProcessor

def runLocal(outf="test.pickle",fileset="test.json",isElectronChannel=True,**options):
    if options['ttbar']:
        p = ttbarProcessor(isElectronChannel,**options)
    else:
        p = MyProcessor(isElectronChannel,**options)

    if options['full']:      
        out = processor.run_uproot_job(
            fileset,
            treename="MuonSystem",
            processor_instance=p,
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
            processor_instance=p,
            executor=processor.iterative_executor,
            executor_args={
                "schema": BaseSchema,
            },
            maxchunks = 1,
            chunksize=100,
        )

    print(out)
    with open(outf,'wb') as f:
        pickle.dump(out,f)
    return

def runLPC(outf="test.pickle",fileset="test.json",isElectronChannel=True,**options):
    import time
    from distributed import Client
    from lpcjobqueue import LPCCondorCluster
    tic = time.time()
    #cluster = LPCCondorCluster(log_directory="/uscms/home/kkwok/log")
    #cluster = LPCCondorCluster()
    #cluster = LPCCondorCluster(shared_temp_directory="/tmp")
    cluster = LPCCondorCluster(shared_temp_directory="/tmp", memory='4GB',
                                 worker_extra_args=['--worker-port 10000:10070', '--nanny-port 10070:10100', '--no-dashboard'],
                                 job_script_prologue=[]) 

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

    if options['ttbar']:
        #p = ttbarProcessor(isElectronChannel,is2017,runSys,saveSkim,debug)
        p = ttbarProcessor(isElectronChannel,**options)
    else:
        p = MyProcessor(isElectronChannel,**options)


    print("Waiting for at least one worker...")
    client.wait_for_workers(4)
    hists, metrics = processor.run_uproot_job(
        fileset,
        treename="MuonSystem",
        processor_instance=p,
        executor=processor.dask_executor,
        executor_args=exe_args,
        chunksize=1000,
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
    parser.add_option('--ttbar', dest='ttbar', action='store_true',default = False, help='Run ttbar proc')
    parser.add_option('--year', dest='year', default = "2018", help='Switch for year dependent pT cut/HLT decision')
    parser.add_option('--full', dest='full', action='store_true',default = False, help='Run full file chunks')
    parser.add_option('--saveSkim', dest='saveSkim', action='store_true',default = False, help='Save skim selections')
    parser.add_option('--debug', dest='debug', action='store_true',default = False, help='run with debug')
    parser.add_option('--runSys', dest='runSys', action='store_true',default = False, help='run systematics variations')
    parser.add_option('--forLimit', dest='forLimit', action='store_true',default = False, help='run only the histograms for limit')
    parser.add_option('--fileset', dest='fileset', default = "test.json", help='input file json')
    parser.add_option('--nJobs', dest='nJobs', default = 4, type=int, help='number of workers in condor')
    parser.add_option('-o', dest='outf', default='HNL_histograms.pickle', help='collection of histograms')
    parser.add_option('-f', dest='inf', default='signals_skim.json' , help='Input fileset')   


    (options, args) = parser.parse_args()
    outf    = options.outf 
    saveSkim    = options.saveSkim
    fileset = options.fileset
    isElectronChannel = not options.muon
    procOptions       = vars(options)
    procOptions       = {k:v for k,v in procOptions.items() if k not in ["fileset","outf","isElectronChannel"]} ## write these 3 option explicitly

    print(" Coffea version = ", coffea.__version__)
    print(" Fileset = ", fileset)
    print(" isElectronChannel = ", isElectronChannel)
    print(" year              = ", options.year)
    print(" runSys            = ", options.runSys)
    print(" outf              = ", outf)
    print(" saveSkim          = ", saveSkim)
    print(" options           = ", procOptions)

    if options.test:
        runLocal("test.pickle","test.json",isElectronChannel,**procOptions)
    elif options.local:
        print("full               = ", options.full)
        runLocal(outf,fileset,isElectronChannel,**procOptions)
    elif options.condor:
        print(" Using nJobs       = ", options.nJobs)
        runLPC(outf,fileset,isElectronChannel,**procOptions)
