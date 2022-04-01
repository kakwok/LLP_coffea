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
    if dataset=="HNL,$m_N$=5" or "HNL" in dataset :
        weights.add("Wpt", compiled["wpt"](gWPt))
    elif dataset=="WJetsToLNu" or  "WJets" in dataset:
        weights.add("Wpt", compiled["wpt_WJ"](gWPt))
        return

def add_ctau_weight(weights,llp_ctau, ctau_old, ctau_new):
    #print("Adding ctau weight")
    w_ctau  =  ak.firsts(np.exp(llp_ctau * (1/ctau_old-1/ctau_new))*(ctau_old/ctau_new))
    if ak.any(np.isnan(w_ctau)):
        w_ctau=ak.nan_to_num(w_ctau,nan=0)  ## zero out nan weights
    weights.add("ctau",w_ctau)
    return

def add_nCluster_weight(weights,name, n_cluster, cls_eff_ratio=1):
    eff = ak.ones_like(n_cluster) * cls_eff_ratio
    #print("Adding ctau weight")
    weights.add(name,np.power(eff,n_cluster))
    return

def load_xsection():
    return compiled['xsections']


def reweightXsec(ctau,mass):
    #ctau in mm
    #mass in GeV
    #xsec in pb
    cof={
        "1p0":9.51561675,
        "2p0":6.04926165,
        "4p0":2.55645182,
        "4p5":1.95648942,
        "5p0":1.42447854,
        "7p0":-0.29378948,
        "10p0":-2.11196473,
    }
    return np.exp(-1*np.log(ctau)+cof[mass])
