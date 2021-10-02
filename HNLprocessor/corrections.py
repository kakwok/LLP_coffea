import numpy as np
import awkward as ak
import gzip
import pickle
import cloudpickle
from coffea.lookup_tools.lookup_base import lookup_base
from coffea import lookup_tools
from coffea import util
import importlib.resources

#with gzip.open("corrections.coffea") as fin:
#    compiled = pickle.load(fin)
#compiled = util.load("/uscms/home/kkwok/work/LLP/CMSSW_10_6_20/src/llp_analyzer/corrections.coffea")
with importlib.resources.path("HNLprocessor","corrections.coffea") as path:
    compiled = util.load(path)

compiled['2018_pileupweight']._values = np.minimum(5, compiled['2018_pileupweight']._values)

def add_pileup_weight(weights, nPU, year='2017', dataset=None):
    if year == '2017' and dataset in compiled['2017_pileupweight_dataset']:
        weights.add(
            'pileup_weight',
            compiled['2017_pileupweight_dataset'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puUp'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puDown'][dataset](nPU),
        )
    else:
        weights.add(
            'pileup_weight',
            compiled[f'{year}_pileupweight'](nPU),
            compiled[f'{year}_pileupweight_puUp'](nPU),
            compiled[f'{year}_pileupweight_puDown'](nPU),
        )

def add_Wpt_kfactor(weights, gWPt, dataset):
    if dataset=="WJetsToLNu" or dataset=="HNL,$m_N$=5" or "HNL" in dataset:
        #print("Adding WpT weight")
        weights.add("Wpt", compiled["wpt"](gWPt))
    else:
        return

def add_ctau_weight(weights,llp_ctau, ctau_old, ctau_new):
    #print("Adding ctau weight")
    weights.add("ctau",ak.firsts(np.exp(llp_ctau * (1/ctau_old-1/ctau_new))*(ctau_old/ctau_new)))
    return

def load_xsection():
    return compiled['xsections']


def reweightXsec(ctau,mass):
    #ctau in mm
    #mass in GeV
    #xsec in pb
    cof={
        "1":9.51561675,
        "2":6.04926165,
        "4":2.55645182,
        "7":-0.29378948,
        "10":-2.11196473,
    }
    return np.exp(-1*np.log(ctau)+cof[mass])
