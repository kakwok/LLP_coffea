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

samples = {}
for era in "A":
    lines  = subprocess.getoutput(
          #'ls /uscms/home/kkwok/eos/llp/SingleMu_2017B/HeavyNeutralLepton_Tree_0.root'
          'ls /uscms/home/kkwok/lpclonglived/HNL/EGamma_2018A/HeavyNeutralLepton_Tree_*.root'
          #'find s /eos/uscms/store/user/lpcbacon/15/JetHTRun2018%s*/*.root -not -empty -ls '%era
          ).split('\n')
    flist = [l.split()[-1] for l in lines]
    if 'directory' in flist:
        flist.remove('directory')
    samples['JetHTRun2017'+era]=flist

treename = "MuonSystem"

#lumivalues = lumi_tools.LumiData("metadata/lumi2016_new.csv")
#lumivalues = lumi_tools.LumiData("metadata/lumi2017.csv.gz")
lumivalues = lumi_tools.LumiData("metadata/lumi2018.csv")
runs  = lumivalues._lumidata[:, 0].astype('u4')
lumis = lumivalues._lumidata[:, 1].astype('u4')
ll = lumi_tools.LumiList(runs,lumis)
print(lumivalues.get_lumi(ll))
#lumivalues = lumi_tools.LumiData("metadata/lumi_2018_golden.csv.gz")

def get_lumilist(dataset, filename):
     print(filename)
     file = uproot.open(filename)
     if treename not in file:
          print("Bad file:", filename)
          return dataset, lumi_tools.LumiList()

     tree = file[treename]
     #run, lumi = tree["run"].array(), tree["luminosityBlock"].array()
     run, lumi = tree["runNum"].array(), tree["lumiSec"].array()
     if len(run)==0: return dataset, lumi_tools.LumiList()
     lumilist = lumi_tools.LumiList(run, lumi)
     return dataset, lumilist


dataset_lumi = {}
nworkers = 12
fileslice = slice(None)
with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = set()
    for dataset, files in samples.items():
        futures.update(executor.submit(get_lumilist, dataset, file) for file in files[fileslice])
    try:
        total = len(futures)
        processed = 0
        while len(futures) > 0:
            finished = set(job for job in futures if job.done())
            for job in finished:
                dataset, accumulator = job.result()
                if dataset in dataset_lumi:
                    dataset_lumi[dataset] += accumulator
                else:
                    dataset_lumi[dataset] = accumulator
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


print("dataset, lumi [/pb], lumisections, unique lumisections")
s = 0
for ds, ll in dataset_lumi.items():
    lumi = lumivalues.get_lumi(ll.array)
    s+=lumi
    nunique = np.unique(ll.array, axis=0).shape[0]
    ntot = ll.array.shape[0]
    print("%50s %0.2f %6d %6d" % (ds, lumi, ntot, nunique))
print("total lumi ",s)

