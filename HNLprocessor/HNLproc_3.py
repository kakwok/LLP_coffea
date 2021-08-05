import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
from coffea import hist, processor
import uproot
# register our candidate behaviors
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema, TreeMakerSchema

from coffea.analysis_tools import Weights, PackedSelection

import HNLprocessor.corrections as corrections 
import warnings
import pickle
import glob

ak.behavior.update(candidate.behavior)


def maskAndFill(denom,selection,value):
    numer = ak.mask(denom,selection)
    numer = ak.fill_none(numer, value) #fill none with same structure
    return ak.flatten(numer)

class MyProcessor(processor.ProcessorABC):
    def __init__(self,debug=False):
        self._debug = debug
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(float),
            "nCluster": hist.Hist("Events",hist.Cat("dataset", "Dataset"),                
                hist.Cat("region", "region"),
                hist.Bin("nCluster", "nCluster", 4, 0, 4),
            ),                   
            "accept": hist.Hist("Events",hist.Cat("dataset", "Dataset"),                
                hist.Bin("gLLP_csc", "gLLP_csc", 2, 0, 2),
            ),                   
            "ClusterSize": hist.Hist("Events",hist.Cat("dataset", "Dataset"),   
                hist.Cat("region", "region"),
                hist.Bin("ClusterSize", "ClusterSize", 50, 0, 2000),
            ),   
            "dphi_cluster_MET": hist.Hist("Events",hist.Cat("dataset", "Dataset"),       
                hist.Cat("region", "region"),                                          
                hist.Bin("ClusterSize", "ClusterSize", 50, 0, 2000),
                hist.Bin("dphi_cluster_MET", r'$\Delta\phi$(cluster,MET)', 50, 0, np.pi),
            ),
            "dphi_cluster_lep": hist.Hist("Events",hist.Cat("dataset", "Dataset"),
                hist.Cat("region", "region"),                                          
                hist.Bin("dphi", r'$\Delta\phi$(cluster,lep)', 50, 0, np.pi),
            ),         
            "metXYCorr": hist.Hist("Events",hist.Cat("dataset", "Dataset"),
                hist.Cat("region", "region"),                                   
                hist.Bin("metXYCorr", "metXYCorr", 50, 0, 500),
            ),

            "MT": hist.Hist("Events",hist.Cat("dataset", "Dataset"),                
                hist.Cat("region", "region"),                            
                hist.Bin("MT", "MT", 50, 0, 200),
            ),            
            "nPU": hist.Hist("Events",hist.Cat("dataset", "Dataset"),                
                hist.Cat("region", "region"),                            
                hist.Bin("nPU", "nPU", 100, 0, 100),
            ),            
            "gWPt": hist.Hist("Events",hist.Cat("dataset", "Dataset"),                
                hist.Cat("region", "region"),                            
                hist.Bin("gWPt", "gWPt", 50,0,500),
            ),            
            "gWPt_noweight": hist.Hist("Events",hist.Cat("dataset", "Dataset"),                
                hist.Cat("region", "region"),                            
                hist.Bin("gWPt", "gWPt", 50,0,500),
            ),            

            "nPU_noweight": hist.Hist("Events",hist.Cat("dataset", "Dataset"),                
                hist.Cat("region", "region"),                            
                hist.Bin("nPU", "nPU", 100, 0, 100),
            ),            

            
        })

    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, events):
        output = self.accumulator.identity()
        dataset = events.metadata['dataset']        
        
        start,stop = events._branchargs['entry_start'],events._branchargs['entry_stop']
        events = uproot.lazy(events._tree)
        events = events[start:stop]
                
        isSignal= ak.any(events.gLLP_csc)        
        isData = not(ak.any(events.gParticleE) or ak.any(events.gLepE))
        output["sumw"][dataset] += len(events)
        if isSignal:        
            csc = ak.any(events.gLLP_csc,axis=1)
            output['accept'].fill(dataset=dataset,gLLP_csc=ak.firsts(events.gLLP_csc)) ## only 1 LLP
            events = events[(csc==1)]        
                                
        llp=ak.zip({
            'pt':events.gLLP_pt,
            'EMfrac':events.gLLP_EMFracE,
            'e':events.gLLP_e,
            'z':events.gLLP_decay_vertex_z ,
            'r':events.gLLP_decay_vertex_r,
        })  
                   
        lep=ak.zip({
            'pt':events.lepPt ,
            'eta':events.lepEta,
            'phi':events.lepPhi,
            'energy':events.lepE,
            'pdgid':events.lepPdgId,
        },with_name='PtEtaPhiELorentzVector',
        behavior=vector.behavior    
        )  
            
        ele   = lep[abs(lep.pdgid)==11]
        muons = lep[abs(lep.pdgid)==13]
        
        cluster= ak.zip(
            {                
                "time":events.cscRechitCluster3Time,
                "size":events.cscRechitCluster3Size,                
                "timeSpread":events.cscRechitCluster3TimeSpread,
                "eta":events.cscRechitCluster3Eta,
                "phi":events.cscRechitCluster3Phi,
                "x":events.cscRechitCluster3X,
                "y":events.cscRechitCluster3Y,
                "z":events.cscRechitCluster3Z,                
                'llp_x':events.cscRechitCluster3_match_gLLP_decay_x,
                'llp_y':events.cscRechitCluster3_match_gLLP_decay_y,
                'llp_z':events.cscRechitCluster3_match_gLLP_decay_z,                 
                "NChamber":events.cscRechitCluster3NChamber,
                "MaxChamber":events.cscRechitCluster3MaxChamber,
                "MaxStation":events.cscRechitCluster3MaxStation,
                "NStation10":events.cscRechitCluster3NStation10,
                "AvgStation10":events.cscRechitCluster3AvgStation10,
                "ME11_12":(events.cscRechitCluster3Me11Ratio+events.cscRechitCluster3Me12Ratio)*events.cscRechitCluster3Size,                
                "llp_match":events.cscRechitCluster3_match_gLLP,
                "RE12":events.cscRechitCluster3_match_RE12_0p4,
                "MB1seg":events.cscRechitCluster3_match_MB1Seg_0p4,
                "RB1":events.cscRechitCluster3_match_RB1_0p4,
                "dphi_cluster_MET":events.cscRechitCluster3MetXYCorr_dPhi,                
            }
        )
        muVeto=ak.zip({
            'e':events.cscRechitCluster3MuonVetoE,
            "time":events.cscRechitCluster3Time,            
            'pt':events.cscRechitCluster3MuonVetoPt,
            'phi':events.cscRechitCluster3MuonVetoPhi,
            'eta':events.cscRechitCluster3MuonVetoEta,
            'LooseIso':events.cscRechitCluster3MuonVetoLooseIso,
            'LooseId':events.cscRechitCluster3MuonVetoLooseId,    
            'TightId':events.cscRechitCluster3MuonVetoTightId,
        })      

        jetVeto=ak.zip({
            'e':events.cscRechitCluster3JetVetoE,
            'pt':events.cscRechitCluster3JetVetoPt,
            'phi':events.cscRechitCluster3JetVetoPhi,
            'eta':events.cscRechitCluster3JetVetoEta,
        })      
        
        ClusterID =((cluster.NStation10>1) & (abs(cluster.eta)<1.9))|\
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==4) & (abs(cluster.eta)<1.8))|\
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==3) & (abs(cluster.eta)<1.6))|\
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==2) & (abs(cluster.eta)<1.6))
        
        muonVeto_mask = ~((muVeto.pt>20) & abs(muVeto.eta<2.4))
        jetVeto_mask  = ~((jetVeto.pt>10)& abs(jetVeto.eta<2.4))
        RE12_veto    = (cluster.RE12==0)
        MB1seg_veto  = (cluster.MB1seg==0)
        RB1_veto     = (cluster.RB1==0)
        ME11_12_veto = (cluster.ME11_12==0)

        OOT_timeCut      = (cluster.time < -12.5) #OOT for data
        IntimeCut      = (cluster.time < 12.5) & (cluster.time>-5) ## In-time otherwise        
        timeSpreadCut= (cluster.timeSpread<20)        

        selection = PackedSelection()        
        
        selection.add('trigger_ele',events.SingleEleTrigger==True)
        selection.add('trigger_mu',events.SingleMuonTrigger==True)
        selection.add('MET',events.metXYCorr>30)
        selection.add('W_CR', (events.MT>70) & (events.MT<90) &(events.metXYCorr>60))
        selection.add('good_electron',ak.num((ele.pt>35) & (abs(ele.eta)<2.4),axis=1)>0)
        
        cls_OOT = cluster[
            ((jetVeto_mask) &(muonVeto_mask))
            & (ME11_12_veto)
            & ((MB1seg_veto) & (RB1_veto))
            & (OOT_timeCut)
            & (timeSpreadCut)
            & (ClusterID)
        ]
        cls_inTime = cluster[
            ((jetVeto_mask) &(muonVeto_mask))
            & (ME11_12_veto)
            & ((MB1seg_veto) & (RB1_veto))
            & (IntimeCut)
            & (timeSpreadCut)
            & (ClusterID)
        ]        
        cls_JetMuVeto = cluster[
            ((jetVeto_mask) &(muonVeto_mask))
        ]        
        cls_JetMuStaVeto = cluster[
            ((jetVeto_mask) &(muonVeto_mask))
            & (ME11_12_veto)
            & ((MB1seg_veto) & (RB1_veto))
        ]        

        selection.add('cls_raw',ak.num(cluster,axis=1)>0)
        selection.add('cls_JetMuVeto',ak.num(cls_JetMuVeto,axis=1)>0)
        selection.add('cls_JetMuStaVeto',ak.num(cls_JetMuStaVeto,axis=1)>0)
        selection.add('cls_InTime',ak.num(cls_inTime,axis=1)>0)
        selection.add('cls_OOT',ak.num(cls_OOT,axis=1)>0)
        selection.add('nJet',events.nJets>0)

        
        regions = {
            "ele_PreSel":['trigger_ele','MET','good_electron'],            
            "ele_W_CR"  :['trigger_ele','MET',"W_CR",'good_electron'],
            "ele_W_CR2"  :['trigger_ele','MET',"W_CR",'good_electron','nJet'],
            "ele_rawcls"  :['trigger_ele','MET','good_electron','cls_raw'],            
            "ele_JetMuVeto"  :['trigger_ele','MET','good_electron','cls_JetMuVeto'],            
            "ele_JetMuStaVeto"  :['trigger_ele','MET','good_electron','cls_JetMuStaVeto'],            
            "ele_SR"    :['trigger_ele','MET','good_electron','cls_InTime'],
            "ele_OOT"   :['trigger_ele','MET','good_electron','cls_OOT'],
            "noselection":[],
        }
        
        ##TODO:  build clusters from passing selections
        cluster_dir= ak.zip(
        {
                'pt':ak.ones_like(events.cscRechitCluster3Eta),
                "eta":events.cscRechitCluster3Eta,
                "phi":events.cscRechitCluster3Phi,
                'mass':ak.zeros_like(events.cscRechitCluster3Eta)
            },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )
        
        cls_lep_pair = ak.cartesian({"cls":cluster_dir,'lep':lep},axis=1,nested=True)
        dphi_lep_cls = cls_lep_pair.cls.delta_phi(cls_lep_pair.lep)        
          
        
            
        for region,cuts in regions.items():
            
            cut = selection.all(*cuts)
            weights = Weights(len(events))
            if not isData:
                corrections.add_pileup_weight(weights, events.npu,'2018')
                corrections.add_Wpt_kfactor(weights, events.gWPt, dataset)
            ##per-event weight
            weight = weights.weight()[cut]
            w_cls      = (weights.weight() * ak.ones_like(cluster.size))[cut] ## use size to pick-up the cluster shape
            w_cls_lep  = (weights.weight() * ak.ones_like(dphi_lep_cls))[cut]

            nCluster = ak.num(cluster[cut],axis=1)
            output["nCluster"].fill(dataset=dataset,region=region,
                                    nCluster=nCluster,
                                   weight=weight)            
            output["ClusterSize"].fill(dataset=dataset,region=region,
                                       ClusterSize=ak.flatten(cluster[cut].size),
                                       weight=ak.flatten(w_cls))        
            output["dphi_cluster_lep"].fill(dataset=dataset,region=region,
                                            dphi =abs(ak.flatten(dphi_lep_cls[cut],axis=None)),
                                            weight=ak.flatten(w_cls_lep,axis=None))
            output["dphi_cluster_MET"].fill(dataset=dataset,region=region,
                                           dphi_cluster_MET=np.abs(ak.flatten(cluster[cut].dphi_cluster_MET)),
                                           ClusterSize=ak.flatten(cluster[cut].size),
                                           weight=ak.flatten(w_cls))

            output["metXYCorr"].fill(dataset=dataset,region=region,
                                     metXYCorr=events[cut].metXYCorr,
                                    weight=weight)                    
            output["MT"].fill(dataset=dataset,region=region,
                              MT=events[cut].MT,
                             weight=weight)        
            output["nPU"].fill(dataset=dataset,region=region,
                              nPU=events[cut].npu,
                             weight=weight)        
            output["gWPt"].fill(dataset=dataset,region=region,
                              gWPt=events[cut].gWPt,
                             weight=weight)        
            output["gWPt_noweight"].fill(dataset=dataset,region=region,
                              gWPt=events[cut].gWPt,
                             )        
            output["nPU_noweight"].fill(dataset=dataset,region=region,
                              nPU=events[cut].npu,
                             )        
        return output

    def postprocess(self, accumulator):
        # set everything to 1/fb scale
        lumi = 1000  # [1/pb]

        scale = {}
        for dataset, dataset_sumw in accumulator['sumw'].items():
            if dataset in corrections.load_xsection():
                scale[dataset] = lumi*corrections.load_xsection()[dataset]/dataset_sumw
            else:
                warnings.warn("Missing cross section for dataset %s.  No normalization applied. " % dataset, RuntimeWarning)
                #scale[dataset] = lumi / dataset_sumw

        for h in accumulator.values():
            if isinstance(h, hist.Hist):
                if self._debug:
                        print("Scaling with scale = " , scale)
                h.scale(scale, axis="dataset")

        return accumulator

