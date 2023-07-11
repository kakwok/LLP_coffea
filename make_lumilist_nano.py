#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import subprocess

import uproot3 as uproot
import numpy as np
from coffea import lumi_tools

def getSamples(sampleInput):
    samples = {}
    #for era in ["A","B","C","D"]:
    for era,cmd in sampleInput.items():
        lines  = subprocess.getoutput(
                cmd
              #'ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleMuon_2018D/HeavyNeutralLepton_Tree_*.root'
              #'find s /eos/uscms/store/user/lpcbacon/15/JetHTRun2018%s*/*.root -not -empty -ls '%era
              ).split('\n')
        flist = [l.split()[-1] for l in lines]
        if 'directory' in flist:
            flist.remove('directory')
        samples[era]=flist
    return samples

def getSamplesFromList(sampleInput):
    samples={}
    for era,txt in sampleInput.items():
        with open(txt,'r') as f:
            flist = f.read().strip().split("\n")
            samples[era] = flist
            print(f'era = {era} , Nfile = ',len(flist))
    return samples

def get_lumilist(dataset, filename,treename="MuonSystem"):
     #print(filename)
     Badfiles = []
     try:
        file = uproot.open(filename)
        if treename not in file:
             print("Bad file:", filename)
             Badfiles.append(filename)
             return dataset, lumi_tools.LumiList(),Badfiles
     except:
        print("Bad file:",filename)
        Badfiles.append(filename)
        return dataset, lumi_tools.LumiList(),Badfiles

     tree = file[treename]
     #run, lumi = tree["run"].array(), tree["luminosityBlock"].array()
     if treename=="MuonSystem":
         run, lumi = tree["runNum"].array(), tree["lumiSec"].array()
     else:
         run, lumi = tree["runNum"].array(), tree["lumiNum"].array()
     if len(run)==0:
         Badfiles.append(filename)
         print("empty file:",filename)
         return dataset, lumi_tools.LumiList(),Badfiles
     lumilist = lumi_tools.LumiList(run, lumi)
     return dataset, lumilist, Badfiles

def printLumi(samples,treename="MuonSystem"):
    if np.all([("2016" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2016_new.csv")
    elif np.all([("2017" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2017.csv.gz")
    elif np.all([("2018" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2018.csv")
    else:
        print(samples.keys()," must contain 2016,2017 or 2018")
        return
    runs  = lumivalues._lumidata[:, 0].astype('u4')
    lumis = lumivalues._lumidata[:, 1].astype('u4')
    ll = lumi_tools.LumiList(runs,lumis)
    print(lumivalues.get_lumi(ll))

    tic = time.time()
    dataset_lumi = {}
    dataset_badFiles = {}
    nworkers = 12
    fileslice = slice(None)
    with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
        futures = set()
        for dataset, files in samples.items():
            futures.update(executor.submit(get_lumilist, dataset, file, treename) for file in files[fileslice])
        try:
            total = len(futures)
            processed = 0
            while len(futures) > 0:
                finished = set(job for job in futures if job.done())
                for job in finished:
                    dataset, accumulator ,badfiles = job.result()
                    if dataset in dataset_lumi:
                        dataset_lumi[dataset] += accumulator
                        dataset_badFiles[dataset] += badfiles
                    else:
                        dataset_lumi[dataset] = accumulator
                        dataset_badFiles[dataset] = badfiles
                    processed += 1
                    if processed % 10 == 0:
                        print("Processing: done with % 4d / % 4d files" % (processed, total))
                futures -= finished
            del finished
        except KeyboardInterrupt:
            print("Ok quitter")
            for job in futures: job.cancel()
        except:
            for job in futures: job.cancel()
            raise
    
    elapsed = time.time() - tic
    print(f"Finished in {elapsed:.1f}s")
   
    fname =[k for k in samples.keys() ][0]
    fout = open("lumi_%s.txt"%fname,"w") 
    print("dataset, lumi [/pb], lumisections, unique lumisections")
    fout.write("dataset, lumi [/pb], lumisections, unique lumisections\n")
    
    s = 0
    for ds, ll in dataset_lumi.items():
        lumi = lumivalues.get_lumi(ll.array)
        s+=lumi
        nunique = np.unique(ll.array, axis=0).shape[0]
        ntot = ll.array.shape[0]
        print("%s %0.2f %6d %6d" % (ds, lumi, ntot, nunique))
        fout.write("%s %0.2f %6d %6d\n" % (ds, lumi, ntot, nunique))
    print("total lumi ",s)
    fout.write("total lumi %s\n"%s)
    print("Bad files:  ")
    for ds, badfiles in dataset_badFiles.items():
        for f in badfiles:
            print(f)
            fout.write("%s\n" % (f))
 
   
if __name__ == '__main__':

    
    EGamma2017 ={
        "EGamma2017B": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleElectron_2017B/HeavyNeutralLepton_Tree_*.root",
        "EGamma2017C": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleElectron_2017C/HeavyNeutralLepton_Tree_*.root",
        "EGamma2017D": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleElectron_2017D/HeavyNeutralLepton_Tree_*.root",
        "EGamma2017E": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleElectron_2017E/HeavyNeutralLepton_Tree_*.root",
        "EGamma2017F": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleElectron_2017F/HeavyNeutralLepton_Tree_*.root",
    }
    #printLumi( getSamples(EGamma2017))

    SingleMuon2017 ={
        "SingleMuon2017B": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleMuon_2017B/HeavyNeutralLepton_Tree_*.root",
        "SingleMuon2017C": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleMuon_2017C/HeavyNeutralLepton_Tree_*.root",
        "SingleMuon2017D": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleMuon_2017D/HeavyNeutralLepton_Tree_*.root",
        "SingleMuon2017E": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleMuon_2017E/HeavyNeutralLepton_Tree_*.root",
        "SingleMuon2017F": "ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleMuon_2017F/HeavyNeutralLepton_Tree_*.root",
    }
    #print(getSamples(SingleMuon2017))
    #printLumi( getSamples(SingleMuon2017))

    Treename = 'ntuples/llp'

    SingleElectron2016 ={
        "SingleElectron_2016B-v1": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016B-v1.txt",
        "SingleElectron_2016B-v2": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016B-v2.txt",
        "SingleElectron_2016C": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016C.txt",
        "SingleElectron_2016D": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016D.txt",
        "SingleElectron_2016E": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016E.txt",
        "SingleElectron_2016F": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016F.txt",
        "SingleElectron_2016G": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016G.txt",
        "SingleElectron_2016H": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016H.txt",
    }
    SingleElectron_2017={
        "SingleElectron_2017B": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017B.txt",
        "SingleElectron_2017C": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017C.txt",
        "SingleElectron_2017D": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017D.txt",
        "SingleElectron_2017E": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017E.txt",
        "SingleElectron_2017F": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017F.txt",
    }
    EGamma_2018 = {
        "EGamma_2018A":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018A.txt",
        "EGamma_2018B":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018B.txt",
        "EGamma_2018C":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018C.txt",
        "EGamma_2018D":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018D.txt",
    }
    #printLumi(getSamplesFromList(SingleElectron2016),Treename)
    #printLumi(getSamplesFromList(SingleElectron_2017),Treename)
    #printLumi(getSamplesFromList(EGamma_2018),Treename)

    SingleMuon_2016={
        "SingleMuon_2016B-v1":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016B-v1.txt",
        "SingleMuon_2016B-v2":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016B-v2.txt",
        "SingleMuon_2016C":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016C.txt",
        "SingleMuon_2016D":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016D.txt",
        "SingleMuon_2016E":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016E.txt",
        "SingleMuon_2016F":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016F.txt",
        "SingleMuon_2016G":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016G.txt",
        "SingleMuon_2016H":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleMuon_2016H.txt",
        }
    SingleMuon_2017={
        "SingleMuon_2017A":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleMuon_2017A.txt",
        "SingleMuon_2017B":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleMuon_2017B.txt",
        "SingleMuon_2017C":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleMuon_2017C.txt",
        "SingleMuon_2017D":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleMuon_2017D.txt",
        "SingleMuon_2017E":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleMuon_2017E.txt",
        "SingleMuon_2017F":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleMuon_2017F.txt",
        }
    SingleMuon_2018={
        "SingleMuon_2018A":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/SingleMuon_2018A.txt",
        "SingleMuon_2018B":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/SingleMuon_2018B.txt",
        "SingleMuon_2018C":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/SingleMuon_2018C.txt",
        "SingleMuon_2018D":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/SingleMuon_2018D.txt",
    }
    #printLumi(getSamplesFromList(SingleMuon_2016),Treename)
    #printLumi(getSamplesFromList(SingleMuon_2017),Treename)
    #printLumi(getSamplesFromList(SingleMuon_2018),Treename)
 
    SingleElectron2016 ={
        "SingleElectron_2016B-v1": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1.txt",
        "SingleElectron_2016B-v2": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2.txt",
        "SingleElectron_2016C"   : "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17.txt",
        "SingleElectron_2016D"   : "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17.txt",
        "SingleElectron_2016E"   : "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17.txt",
        "SingleElectron_2016F"   : "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17.txt",
        "SingleElectron_2016G"   : "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17.txt",
        "SingleElectron_2016H"   : "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17.txt",
    }
    SingleElectron_2017={
        "SingleElectron_2017B": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017.txt",
        "SingleElectron_2017C": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017.txt",
        "SingleElectron_2017D": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017.txt",
        "SingleElectron_2017E": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017.txt",
        "SingleElectron_2017F": "../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleElectron/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017F-17Nov2017.txt",
    }
    EGamma_2018 = {
        "EGamma_2018A"      :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/EGamma/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v2.txt",
        "EGamma_2018A_part2":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/EGamma/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018-v2_part2.txt",
        "EGamma_2018B"      :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/EGamma/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018.txt",
        "EGamma_2018C"      :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/EGamma/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018.txt",
        "EGamma_2018D"      :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D_AOD/v5/sixie/EGamma/Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco.txt",
    }
    #printLumi(getSamplesFromList(SingleElectron2016),Treename)
    #printLumi(getSamplesFromList(SingleElectron_2017),Treename)
    printLumi(getSamplesFromList(EGamma_2018),Treename)

    print("Lumi for single muon")
    SingleMuon_2016={
        "SingleMuon_2016B-v1":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver1.txt",
        "SingleMuon_2016B-v2":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016B-07Aug17_ver2.txt",
        "SingleMuon_2016C"   :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016C-07Aug17.txt",
        "SingleMuon_2016D"   :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016D-07Aug17.txt",
        "SingleMuon_2016E"   :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016E-07Aug17.txt",
        "SingleMuon_2016F"   :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016F-07Aug17.txt",
        "SingleMuon_2016G"   :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016G-07Aug17.txt",
        "SingleMuon_2016H"   :"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2016_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2016_AOD_Run2016H-07Aug17.txt",
        }
    SingleMuon_2017={
        "SingleMuon_2017B":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017B-17Nov2017.txt",
        "SingleMuon_2017C":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017C-17Nov2017.txt",
        "SingleMuon_2017D":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017D-17Nov2017.txt",
        "SingleMuon_2017E":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017E-17Nov2017.txt",
        "SingleMuon_2017F":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2017_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2017_AOD_Run2017F-17Nov2017.txt",
        }
    SingleMuon_2018={
        "SingleMuon_2018A":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018A-17Sep2018.txt",
        "SingleMuon_2018B":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018B-17Sep2018.txt",
        "SingleMuon_2018C":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018ABC_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2018_17Sept2018_AOD_Run2018C-17Sep2018.txt",
        "SingleMuon_2018D":"../llp_analyzer/lists/displacedJetMuonNtuple/V1p17/Data2018D_AOD/v5/sixie/SingleMuon/Run2_displacedJetMuonNtupler_V1p17_Data2018D_17Sept2018_AOD_Run2018D-PromptReco.txt",
    }
    printLumi(getSamplesFromList(SingleMuon_2016),Treename)
    printLumi(getSamplesFromList(SingleMuon_2017),Treename)
    printLumi(getSamplesFromList(SingleMuon_2018),Treename)

  

# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016B-v1.txt |  390 ++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016B-v2.txt | 3212 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016C.txt    | 1035 +++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016D.txt    | 1551 +++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016E.txt    | 1384 +++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016F.txt    |  994 ++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016G.txt    | 2348 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2016/MuonHitsOnly/SingleElectron_2016H.txt    | 2727 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017A.txt    |  126 +++
# lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017B.txt    | 1396 +++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017C.txt    | 2985 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017D.txt    | 1448 ++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017E.txt    | 2301 ++++++++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2017/MuonHitsOnly/SingleElectron_2017F.txt    | 3060 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018A.txt            | 2060 +++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018B.txt            | 2015 ++++++++++++++++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018C.txt            | 1428 ++++++++++++++++++++++++++++++++++
# lists/displacedJetMuonNtuple/V1p17/Data2018/MuonHitsOnly/EGamma_2018D.txt 


