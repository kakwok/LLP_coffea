import numpy as np
import awkward as ak
import gzip
import pickle
import cloudpickle
from coffea.lookup_tools.lookup_base import lookup_base
from coffea import lookup_tools
from coffea import util
import importlib.resources
from coffea.lookup_tools import extractor
import HNLprocessor.util as HNLutil
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
        weights.add("Wpt", compiled["wpt"](gWPt)
                            , compiled["wptUp"](gWPt)
                            , compiled["wptDown"](gWPt))
    elif dataset=="WJetsToLNu" or  "WJets" in dataset:
        weights.add("Wpt"    , compiled["wpt_WJ"](gWPt)
                             ,compiled["wpt_WJUp"](gWPt)
                             ,compiled["wpt_WJDown"](gWPt))
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


def reweightXsec(ctau,mass,isTau=False):
    #ctau in mm
    #mass in GeV
    #xsec in pb
    if isTau:
        return  HNLutil.f_xsec_tau(mass)(ctau)
    else:
        return  HNLutil.f_xsec(mass)(ctau)

# Root files from https://twiki.cern.ch/twiki/bin/view/CMS/MuonLegacy2018#Medium_pT_from_15_to_120_GeV
def add_muonSFs(weights, leadingmuon):

    for sf in compiled['muonsf_keys']:

        lep_pt = np.array(ak.fill_none(leadingmuon.pt, 0.))
        lep_eta = np.array(ak.fill_none(leadingmuon.eta, 0.))

        if 'value' in sf:
            
            nom = compiled['muonsf_evaluator'][sf](np.abs(lep_eta),lep_pt)
            shift = compiled['muonsf_evaluator'][sf.replace('_value','_error')](np.abs(lep_eta),lep_pt)
                      
            weights.add(sf, nom, shift, shift=True)

# Root files from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaRunIIRecommendations#Electron_Scale_Factors
def add_electronSFs(weights, leadingelectron,year):

    lep_pt = np.array(ak.fill_none(leadingelectron.pt, 0.))
    lep_eta = np.array(ak.fill_none(leadingelectron.eta, 0.))

    nom = compiled['elesf_evaluator']["electron_SF_2018_value"](np.abs(lep_eta),lep_pt)
    shift = compiled['elesf_evaluator']["electron_SF_2018_error"](np.abs(lep_eta),lep_pt)
    weights.add("electron_SF_2018_value", nom, shift, shift=True)
    ## add trigger SF
    nom   = compiled['elesf_evaluator']["electron_trigger_SF_%s_value"%year](np.abs(lep_eta),lep_pt)
    shift = compiled['elesf_evaluator']["electron_trigger_SF_%s_error"%year](np.abs(lep_eta),lep_pt)
    
    weights.add("electron_trigger_SF_value", nom, shift, shift=True)
    return
    

def build_lumimask(filename):
    from coffea.lumi_tools import LumiMask
    with importlib.resources.path("HNLprocessor.data", filename) as path:
        return LumiMask(path)


lumiMasks = {
    "2016": build_lumimask("Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
    "2017": build_lumimask("Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
    "2018": build_lumimask("Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
}
