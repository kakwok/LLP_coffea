import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
from coffea import hist, processor
import uproot
# register our candidate behaviors
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema, TreeMakerSchema
from coffea.lumi_tools import LumiMask
import time

from coffea.analysis_tools import Weights, PackedSelection

import HNLprocessor.corrections as corrections 
from  HNLprocessor.corrections import lumiMasks
from  HNLprocessor.histograms import histograms 
from  HNLprocessor import util 
import warnings
import pickle
import glob
import os
import XRootD
import XRootD.client


ak.behavior.update(candidate.behavior)

def uproot_writeable(events, keys):
    ev = {}
    for bname in keys:
        if not events[bname].fields:
            ev[bname] = ak.packed(ak.without_parameters(events[bname]))
        else:
            ev[bname] = ak.zip(
                {
                    n: ak.packed(ak.without_parameters(events[bname][n]))
                    for n in events[bname].fields
                    if is_rootcompat(events[bname][n])
                }
            )
    return ev


def maskAndFill(denom,selection,value):
    numer = ak.mask(denom,selection)
    numer = ak.fill_none(numer, value) #fill none with same structure
    return ak.flatten(numer)

# e.g. pair_dict = {"cls":cluster,"bjets":bjets}
def delta_r_pairs(pair_dict):
    pairs = ak.cartesian(pair_dict,axis=1,nested=True)
    main_obj = getattr(pairs,pairs.fields[0])
    sec_obj = getattr(pairs,pairs.fields[1])    
    dr_pairs = main_obj.delta_r(sec_obj)       
    return dr_pairs,pairs
def delta_phi_pairs(pair_dict):
    pairs = ak.cartesian(pair_dict,axis=1,nested=True)
    main_obj = getattr(pairs,pairs.fields[0])
    sec_obj = getattr(pairs,pairs.fields[1])    
    dphi_pairs = main_obj.delta_phi(sec_obj)       
    return dphi_pairs,pairs

def buildMask(allMasks,cutnames):
    if type(allMasks)==type({}):
        ## build masks with event allMasks[cutnames]
       allcuts = allMasks[cutnames[0]]
       for i,cutname in enumerate(cutnames):
           allcuts = (allcuts) & allMasks[cutname]
    else:
        ## build masks with event zip.cutnames
        allcuts = getattr(allMasks,cutnames[0])
        for s in cutnames:
            allcuts = (allcuts) & getattr(allMasks,s)
    return allcuts

class MyProcessor(processor.ProcessorABC):
    """
    Processor for HNL analysis with muon shower in the muon system. Make histograms from custom flat input tree.

    Parameters 
    ----------
    isElectronChannel : bool (default is True)
        If true, make electron channel selections; else make muon channel selections;
    year : str (default is "2018")
        Switches for trigger decisions and prompt lepton pT
    runSys : bool (default is False)
        Run systematic variations of signals
    forLimit : bool (default is False)
        Produce only the histograms that are used for limit setting
    saveSkim : bool (default is False)
        Output events passing selections with reduced set of info. Modify processor to configure selection.
    debug   : bool (default is False)
        Print debug information.
    
    """
    def __init__(self,isElectronChannel=True,**options):
        defaultOptions = { 'debug': False, 'saveSkim': False, 'runSys':False,"year":"2018","forLimit":False}
        options = { **defaultOptions, **options }
        self._debug = options['debug']
        self._saveSkim = options['saveSkim']
        self.isElectronChannel = isElectronChannel
        self.isMuonChannel = not(isElectronChannel)
        self.llp = None
        self._year = options['year']
        self._runSys = options['runSys']
        self._forLimit = options['forLimit']
        ##define histograms 
        histograms['sumw']= processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator( histograms )

    @property
    def accumulator(self):
        return self._accumulator

    def buildRecoMuons(self,events):
        muons = util.pack(events,"muon")
        muons = ak.with_field(muons,muons.Phi,"phi")
        muons = ak.with_field(muons,muons.Eta,"eta")
        muons = muons[ak.argsort(muons.Pt, axis=1,ascending=False)]

        #glb_muons = muons[muons.IsGlobal==True]
        #glb_muons_pt2 = glb_muons[glb_muons.Pt>2]
        #all_muons_pt5 =muons[muons.Pt>5.0]
        
        #return muons#,glb_muons_pt2,all_muons_pt5
        return muons

    def buildGoodLeptons(self,events):
        lep=ak.zip({
            'pt':events.lepPt ,
            'eta':events.lepEta,
            'phi':events.lepPhi,
            'energy':events.lepE,
            'pdgid':events.lepPdgId,
            'passId':events.lepPassId,
        },with_name='PtEtaPhiELorentzVector',
        behavior=vector.behavior    
        )  
            
        ele   = lep[abs(lep.pdgid)==11]
        muons = lep[abs(lep.pdgid)==13]
        if self._year=="2016":
            good_ele = ele[(ele.pt>30) & (abs(ele.eta)<2.5) & (ele.passId)]
        elif self._year =="2017":
            good_ele = ele[(ele.pt>35) & (abs(ele.eta)<2.5) & (ele.passId)]
        elif self._year =="2018":
            #good_ele = ele[(ele.pt>32) & (abs(ele.eta)<2.5) & (ele.passId)]
            good_ele = ele[(ele.pt>35) & (abs(ele.eta)<2.5) & (ele.passId)]
            
        if self._year=="2017":
            good_mu  = muons[(muons.pt>28)&(abs(muons.eta)<2.4) & (muons.passId)]   ## use pT>28 GeV for 2017
        else:
            good_mu  = muons[(muons.pt>25)&(abs(muons.eta)<2.4) & (muons.passId)] 
 
        if self.isElectronChannel:  good_lep = good_ele
        elif self.isMuonChannel: good_lep = good_mu
        return good_lep,ele, muons

    def buildLLP(self,events):
        llp=ak.zip({
            'pt':events.gLLP_pt,
            #'EMfrac':events.gLLP_EMFracE,
            'e':events.gLLP_e,
            'eta':events.gLLP_eta,
            'phi':events.gLLP_phi,
            'z':events.gLLP_decay_vertex_z ,
            'r':events.gLLP_decay_vertex_r,
            'ctau':events.gLLP_ctau,
        })
        llp['csc_loose']= (abs(llp.eta)<2.4) &(abs(llp.eta)>1.9)& (llp.r<695.5) &( abs(llp.z)>400) &( abs(llp.z)<1100)
        llp['dt_loose'] = (abs(llp.z)<661) & (llp.r>200) & (llp.r<800) 
        llp['csc_tight']= (abs(llp.eta)<2.4) &(abs(llp.eta)>1.9)& (llp.r<695.5) &( abs(llp.z)>730) &( abs(llp.z)<1100)
        llp['dt_tight'] = (abs(llp.z)<661) & (llp.r>460) & (llp.r<800) 
        self.llp = llp
        return llp 

    def buildGenParticles(self,events):
        gParticle = ak.zip({
            "E":events.gParticleE,
            "Eta":events.gParticleEta,
            "Phi":events.gParticlePhi,    
            "Id":events.gParticleId,
            "MotherId":events.gParticleMotherId,
            "MotherIndex":events.gParticleMotherIndex,
            "Pt":events.gParticlePt,
            "Status":events.gParticleStatus,
            "ProdVertexX":events.gParticleProdVertexX,
            "ProdVertexY":events.gParticleProdVertexY,
            "ProdVertexZ":events.gParticleProdVertexZ,
        })
        return gParticle
        
    def buildCSCcluster(self, events,good_lep):
        cluster_dir= ak.zip(
        {
                'pt':ak.ones_like(events.cscRechitCluster3Eta),
                "eta":events.cscRechitCluster3Eta,
                "phi":events.cscRechitCluster3Phi,
                'mass':ak.zeros_like(events.cscRechitCluster3Eta)
            },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )
        #compute dphi with selected electron with highest pT
        dphi_cluster_lep = ak.fill_none(cluster_dir.delta_phi(ak.firsts(good_lep)),-999,axis=None)
        dr_cluster_lep = ak.fill_none(cluster_dir.delta_r(ak.firsts(good_lep)),-999,axis=None)
        #dphi_cluster_mu = ak.fill_none(cluster_dir.delta_phi(ak.firsts(muons)),-999)
      
        cluster= ak.zip(
            {                
                "time":events.cscRechitCluster3TimeTotal,
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
                'MuonVetoPt':events.cscRechitCluster3MuonVetoPt,
                'MuonVetoEta':events.cscRechitCluster3MuonVetoEta,
                'JetVetoPt':events.cscRechitCluster3JetVetoPt,
                'JetVetoEta':events.cscRechitCluster3JetVetoEta,
                "dphi_cluster_MET":events.cscRechitCluster3MetXYCorr_dPhi,                
                "dphi_cluster_lep":dphi_cluster_lep,                
                "dr_cluster_lep":dr_cluster_lep,                
            },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )
        #dr_glbMuPt2 = ak.fill_none(cluster.delta_r(ak.firsts(glb_muons_pt2)),1000,axis=None)
        #dr_allMuPt5 = ak.fill_none(cluster.delta_r(ak.firsts(all_muons_pt5)),1000,axis=None)
        
        #cluster = ak.with_field(cluster,dr_glbMuPt2,"dr_glbMuPt2")
        #cluster = ak.with_field(cluster,dr_allMuPt5,"dr_allMuPt5")

        return cluster

    def selectCSCcluster(self,cluster,events):
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
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==4) & (abs(cluster.eta)<1.6))|\
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==3) & (abs(cluster.eta)<1.6))|\
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==2) & (abs(cluster.eta)<1.6))
        
        llp_match = (cluster.llp_match==True) 
        muonVeto = ~((muVeto.pt>20) & abs(muVeto.eta<2.4))
        jetVeto  = ~((jetVeto.pt>10)& abs(jetVeto.eta<3.0))
        RE12_veto     = (cluster.RE12==0)
        MB1seg_veto   = (cluster.MB1seg==0)
        RB1_veto      = (cluster.RB1==0)
        ME11_12_veto  = (cluster.ME11_12==0)

        oot_timecut   = (cluster.time < -12.5) #OOT for data
        IntimeCut     = (cluster.time < 12.5) & (cluster.time>-5) ## In-time otherwise        
        timeSpreadCut = (cluster.timeSpread<20)        
        dphi_met      = (abs(cluster.dphi_cluster_MET)<0.7)        
        dphi_lep      = (abs(cluster.dphi_cluster_lep)>2.5)      
        dr_lep      = (cluster.dr_cluster_lep>0.8)
        #dr_glbMu      = (cluster.dr_glbMuPt2>0.8)
        #dr_allMu      = (cluster.dr_allMuPt5>0.8)

        clusterMasks = ak.zip({
            "ClusterID"     : ClusterID     ,  
            "muonVeto" : muonVeto ,
            "jetVeto"  : jetVeto  ,
            "RE12_veto"     : RE12_veto     ,
            "MB1seg_veto"   : MB1seg_veto   ,
            "RB1_veto"      : RB1_veto      ,
            "ME11_12_veto"  : ME11_12_veto  ,
                                         
            "OOT_timeCut"   : oot_timecut   ,
            "IntimeCut"     : IntimeCut     ,
            "timeSpreadCut" : timeSpreadCut ,
            "dphi_MET"      : dphi_met      ,
            "dphi_lep"      : dphi_lep      ,
            "dr_lep"      : dr_lep      ,
            "llp_match"      : llp_match      ,
            #"dr_glbMu"      : dr_glbMu      ,
            #"dr_allMu"      : dr_allMu      ,
            })
        return clusterMasks 


    def buildDTcluster(self, events,good_lep):
        dt_cluster_dir= ak.zip(
        {
                'pt':ak.ones_like(events.dtRechitClusterEta),
                "eta":events.dtRechitClusterEta,
                "phi":events.dtRechitClusterPhi,
                'mass':ak.zeros_like(events.dtRechitClusterEta)
            },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )
        dphi_dt_cluster_lep = ak.fill_none(dt_cluster_dir.delta_phi(ak.firsts(good_lep)),-999,axis=None)
        dr_dt_cluster_lep = ak.fill_none(dt_cluster_dir.delta_r(ak.firsts(good_lep)),-999,axis=None)
        dt_cluster = ak.zip(
            {
                 "size":events.dtRechitClusterSize,
                 "X":events.dtRechitClusterX,
                 "Y":events.dtRechitClusterY,
                 "Z":events.dtRechitClusterZ,
                 "R":np.sqrt(events.dtRechitClusterX**2+events.dtRechitClusterY**2),
                 "llp_x":events.dtRechitCluster_match_gLLP_decay_x,
                 "llp_y":events.dtRechitCluster_match_gLLP_decay_y,
                 "llp_z":events.dtRechitCluster_match_gLLP_decay_z,
                 "llp_match":events.dtRechitCluster_match_gLLP,
                 "eta":events.dtRechitClusterEta,
                 "phi":events.dtRechitClusterPhi,
                 "JetVetoPt":events.dtRechitClusterJetVetoPt,
                 "JetVetoEta":events.dtRechitClusterJetVetoEta,
                 "MuonVetoPt":events.dtRechitClusterMuonVetoPt,
                 "MuonVetoLooseId":events.dtRechitClusterMuonVetoLooseId,
                 "MuonVetoGlobal":events.dtRechitClusterMuonVetoGlobal,
                 "NStation10":events.dtRechitClusterNStation10,
                 "AvgStation10":events.dtRechitClusterAvgStation10,
                 "MaxStation":events.dtRechitClusterMaxStation,
                 "nRPC":events.dtRechitCluster_match_RPChits_dPhi0p5,
                 "nMB1_cosmic_minus":events.dtRechitCluster_match_MB1hits_cosmics_minus,
                 "nMB1_cosmic_plus":events.dtRechitCluster_match_MB1hits_cosmics_plus,
                 "nMB1":events.dtRechitCluster_match_MB1hits_0p5,
                 "nMB2":events.dtRechitClusterNSegmentStation2,
                 "nMB3":events.dtRechitClusterNSegmentStation3,
                 "nMB4":events.dtRechitClusterNSegmentStation4,
                 "rpcBx":events.dtRechitCluster_match_RPCBx_dPhi0p5,
                 "dphi_cluster_MET":events.dtRechitClusterMetEENoise_dPhi,
                 "dphi_cluster_lep":dphi_dt_cluster_lep,
                 "dr_cluster_lep":dr_dt_cluster_lep,
            },with_name="PtEtaPhiMLorentzVector",
             behavior=vector.behavior,
        )
        #dr_glbMuPt2 = ak.fill_none(dt_cluster.delta_r(ak.firsts(glb_muons_pt2)),1000,axis=None)# fill 1000 if no glb muons
        #dr_allMuPt5 = ak.fill_none(dt_cluster.delta_r(ak.firsts(all_muons_pt5)),1000,axis=None)# fill 1000 if no  muons

        #dt_cluster = ak.with_field(dt_cluster,dr_glbMuPt2,"dr_glbMuPt2")
        #dt_cluster = ak.with_field(dt_cluster,dr_allMuPt5,"dr_allMuPt5")

        eta_0 = ak.full_like(events.weight,0.3,dtype=float)
        eta_1 = ak.full_like(events.weight,-0.3,dtype=float)

        phi_0 = ak.full_like(events.weight,1.7,dtype=float)
        phi_1 = ak.full_like(events.weight,1.15,dtype=float)

        deadzone_1 = ak.zip(
            {
            'pt':ak.ones_like(events.weight),
            "eta":eta_0,
            "phi":phi_0,
            'mass':ak.ones_like(events.weight)
            },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )

        deadzone_2 = ak.zip(
        {
            'pt':ak.ones_like(events.weight),
            "eta":eta_1,
            "phi":phi_1,
            'mass':ak.ones_like(events.weight)
            },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )

        dr_dt_cluster_dz1 = ak.fill_none(dt_cluster.delta_r(deadzone_1),-999,axis=None)

        dt_cluster = ak.with_field(dt_cluster,dr_dt_cluster_dz1<0.4,"Deadzone_1")

        dr_dt_cluster_dz2 = ak.fill_none(dt_cluster.delta_r(deadzone_2),-999,axis=None)

        dt_cluster = ak.with_field(dt_cluster,dr_dt_cluster_dz2<0.4,"Deadzone_2")

        dt_IsNoise = (events.dtRechitClusterMaxStation==2)&(events.dtRechitClusterWheel==1)\
                            &(events.dtRechitClusterPhi<=0.7)&(events.dtRechitClusterPhi>=0.4)\
                            &( ((events.runNum>=275750)&(events.runNum<=275950))|\
                                ((events.runNum>=274150)&(events.runNum<=274200)) )  #2016C,2016B
    
        dt_cluster = ak.with_field(dt_cluster,dt_IsNoise==True,"IsNoise")

        return dt_cluster
 
    def selectDTcluster(self,dt_cluster,events):
        dt_jetVeto  = ~((dt_cluster.JetVetoPt>20.0) & (abs(dt_cluster.JetVetoEta)<3.0))
        #dt_muonVeto = ~( (dt_cluster.MuonVetoPt>10.0) & (dt_cluster.MuonVetoLooseId==True))
        dt_muonVeto = ~( (dt_cluster.MuonVetoPt>10.0) )
        dt_MB1veto  = (dt_cluster.nMB1<=1)
        dt_etaveto  = ~((abs(dt_cluster.eta)<=0.9) & (abs(dt_cluster.eta)>=0.78))
        dt_RPC      = (dt_cluster.nRPC>=1)
        dt_MB1adj   = (dt_cluster.nMB1_cosmic_minus<=8) & (dt_cluster.nMB1_cosmic_plus<=8)
        dt_time     = (dt_cluster.rpcBx==0)
        dt_OOT      = (dt_cluster.rpcBx>=-100)&(dt_cluster.rpcBx<0)
        dt_dphi_MET  = (abs(dt_cluster.dphi_cluster_MET)<0.7)
        dt_size      = (dt_cluster.size>=100)
        dr_lep      = (dt_cluster.dr_cluster_lep>0.8)
        #dr_glbMu      = (dt_cluster.dr_glbMuPt2>0.8)
        #dr_allMu      = (dt_cluster.dr_allMuPt5>0.8)
        dt_deadzones = ~(dt_cluster.Deadzone_1) & ~(dt_cluster.Deadzone_2)
        dt_noise      = (dt_cluster.IsNoise==False) 
        dt_maxMB2      = (dt_cluster.MaxStation==2) 
        dt_maxMB34     = (dt_cluster.MaxStation>=3) 
        dt_llp_match     = (dt_cluster.llp_match==True) 
        clusterMasks = ak.zip({
                "dt_jetVeto"  :dt_jetVeto  ,
                "dt_muonVeto" :dt_muonVeto ,
                "dt_MB1veto"  :dt_MB1veto  ,
                "dt_etaveto"  :dt_etaveto  ,
                "dt_RPC"      :dt_RPC      ,
                "dt_MB1adj"   :dt_MB1adj   ,
                "dt_time"     :dt_time     ,
                "dt_OOT"      :dt_OOT      ,
                "dt_dphi_MET" :dt_dphi_MET ,
                "dt_size"     :dt_size     ,
                "dr_lep"      :dr_lep     ,
                #"dr_glbMu"     :dr_glbMu     ,
                #"dr_allMu"     :dr_allMu     ,
                "dt_deadzones": dt_deadzones,
                "dt_maxMB2"   : dt_maxMB2,
                "dt_maxMB34"  : dt_maxMB34,
                "dt_llp_match"  : dt_llp_match,
                "dt_noise"    : dt_noise,
        })
        return clusterMasks
    
    def buildHLTdecision(self,events):
        hlt = events.HLTDecision    
        ## from trigger_names_llp_v3
        bitmap = {
             "HLT_Ele27_WPTight_Gsf"                  :   79,  
             "HLT_Ele32_WPTight_Gsf"                  :   87, 
             "HLT_Ele35_WPTight_Gsf"                  :  625,
             "HLT_Ele115_CaloIdVT_GsfTrkIdT"          :  521,
             "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"  :  516,
             "HLT_Photon175"                          :  409,
             "HLT_Photon200"                          :  775,
             ## Muon triggers                                                      
             "HLT_IsoMu24"                            :  135,
             "HLT_IsoTkMu24"                          :  141,
             "HLT_IsoMu27"                            :  136,
             "HLT_Mu50"                               :  196,
             "HLT_OldMu100"                           :  663,
             "HLT_TkMu100"                            :  664,
             "HLT_TkMu50"                             :  197,
        }
        hlt_result = ak.Array([False]*len(hlt))
        if self.isElectronChannel:
            if self._year=="2016":
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Ele27_WPTight_Gsf"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Ele115_CaloIdVT_GsfTrkIdT"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Photon175"]])
            elif self._year=="2017":
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Ele35_WPTight_Gsf"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Ele115_CaloIdVT_GsfTrkIdT"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Photon200"]])
            elif self._year=="2018":
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Ele32_WPTight_Gsf"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Ele115_CaloIdVT_GsfTrkIdT"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Photon200"]])
            else:
                warnings.warn(" Invalid year is selected. No HLT decision is selected"  )
                return np.ones(len(events))
        else:
            if self._year=="2016":
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_IsoMu24"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_IsoTkMu24"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Mu50"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_TkMu50"]])
            elif self._year=="2017":
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_IsoMu27" ]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Mu50"    ]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_OldMu100"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_TkMu100" ]])
            elif self._year=="2018":
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_IsoMu24" ]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_Mu50"    ]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_OldMu100"]])
                hlt_result = (hlt_result) | (hlt[:,bitmap["HLT_TkMu100" ]])
            else:
                warnings.warn(" Invalid year is selected. No HLT decision is selected"  )
                return np.ones(len(events))
        if self._debug:
            print(" HLT pass = ", ak.sum(hlt_result))
        return hlt_result

        
    def buildSelectionMasks(self,events,good_lep,cluster,clusterMasks,dt_cluster,dt_clusterMasks):
        selectionMasks =   {}

        selectionMasks['Acceptance_csc']   =ak.firsts(events.gLLP_csc)==1
        selectionMasks['Acceptance_dt']=ak.firsts(events.gLLP_dt)==1
        selectionMasks['Acceptance_csc_loose']=ak.firsts(self.llp.csc_loose)==1
        selectionMasks['Acceptance_dt_loose']=ak.firsts(self.llp.dt_loose)==1
        selectionMasks['Acceptance_csc_tight']=ak.firsts(self.llp.csc_tight)==1
        selectionMasks['Acceptance_dt_tight']=ak.firsts(self.llp.dt_tight)==1
        selectionMasks['METfilters']   =events.Flag2_all==True
        #selectionMasks['trigger_ele']  =events.SingleEleTrigger==True
        #selectionMasks['trigger_mu']   =events.SingleMuonTrigger==True
        selectionMasks['trigger_ele']  = self.buildHLTdecision(events)
        selectionMasks['trigger_mu']   = self.buildHLTdecision(events)
        selectionMasks['good_lepton']  =ak.num(good_lep,axis=1)==1
        selectionMasks['MET']          =events.metEENoise>=30
        selectionMasks['n_cls']        =ak.num(cluster,axis=1)>=1
        selectionMasks['n_cls_dt']     =ak.num(dt_cluster,axis=1)>=1

        clusterMasks["neg_ME11_12_veto"] = ~clusterMasks['ME11_12_veto']  #make veto mask
        CSC_sel_ABCD = ["ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto","RE12_veto","IntimeCut","timeSpreadCut","ClusterID"]
        CSC_sel_OOT  = ["ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto","RE12_veto",
                        "OOT_timeCut","timeSpreadCut","ClusterID"]
        CSC_sel_negME11 = ["neg_ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto", "IntimeCut","timeSpreadCut","ClusterID"]

        selectionMasks['cls_PreSel']  = buildMask(clusterMasks,["llp_match"])         
        selectionMasks['cls_ABCD']  = buildMask(clusterMasks,CSC_sel_ABCD)
        selectionMasks['cls_OOT']   = buildMask(clusterMasks,CSC_sel_OOT)
        selectionMasks['cls_negME11']   = buildMask(clusterMasks,CSC_sel_negME11)

        dt_clusterMasks["neg_dt_MB1veto"] = ~dt_clusterMasks["dt_MB1veto"] #make veto dt mask
        selectionMasks['cls_StatVeto']     =  buildMask(clusterMasks,['ME11_12_veto','MB1seg_veto','RB1_veto',"RE12_veto"])     
        selectionMasks['cls_JetMuVeto']    =  buildMask(clusterMasks,['jetVeto','muonVeto'])                
        selectionMasks['cls_JetMuStaVeto'] =  buildMask(clusterMasks,['jetVeto','muonVeto','ME11_12_veto','MB1seg_veto','RB1_veto',"RE12_veto"])
        selectionMasks['MCclosure'] =  buildMask(clusterMasks,['dr_lep','jetVeto','muonVeto',"IntimeCut"])

        DT_sel_OOT      = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_OOT" ]
        DT_sel_OOT_MB2      = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_maxMB2","dt_OOT" ]
        DT_sel_OOT_MB34      = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_maxMB34","dt_OOT" ]
        DT_sel_ABCD     = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_time"]
        DT_sel_ABCD_MB2 = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_maxMB2","dt_time"]
        DT_sel_ABCD_MB34 = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_maxMB34","dt_time"]

        #DT_sel_ABCD     = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj",'dt_etaveto',"dt_deadzones","dt_noise","dt_time"]
        #DT_sel_ABCD_MB2 = ["dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj",'dt_etaveto',"dt_deadzones","dt_noise","dt_maxMB2","dt_time"]
        #DT_sel_ABCD_MB34 = ["dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj",'dt_etaveto',"dt_deadzones","dt_noise","dt_maxMB34","dt_time"]

        DT_sel_vetos = ["dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_deadzones","dt_noise"          ]
        DT_sel_negMB1      = ["neg_dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_time"]
        DT_sel_negMB1_MB2  = ["neg_dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_maxMB2","dt_time"]
        DT_sel_negMB1_MB34 = ["neg_dt_MB1veto","dt_jetVeto","dt_muonVeto" ,"dt_RPC","dt_MB1adj","dt_deadzones","dt_noise","dt_maxMB34","dt_time"]
        

        selectionMasks['dt_cls_PreSel']  = buildMask(dt_clusterMasks,["dt_llp_match"])         
        selectionMasks['dt_cls_OOT']  = buildMask(dt_clusterMasks,DT_sel_OOT)         
        selectionMasks['dt_cls_OOT_MB2']  = buildMask(dt_clusterMasks,DT_sel_OOT_MB2)         
        selectionMasks['dt_cls_OOT_MB34']  = buildMask(dt_clusterMasks,DT_sel_OOT_MB34)         
        selectionMasks['dt_cls_ABCD']  = buildMask(dt_clusterMasks,DT_sel_ABCD)         
        selectionMasks['dt_cls_ABCD_MB2'] = buildMask(dt_clusterMasks,DT_sel_ABCD_MB2)         
        selectionMasks['dt_cls_ABCD_MB34'] = buildMask(dt_clusterMasks,DT_sel_ABCD_MB34)         
        selectionMasks['dt_JetMuStaVeto'] =  buildMask(dt_clusterMasks,DT_sel_vetos)
        selectionMasks['dt_cls_negMB1'] = buildMask(dt_clusterMasks,DT_sel_negMB1)         
        selectionMasks['dt_cls_negMB1_MB2'] = buildMask(dt_clusterMasks,DT_sel_negMB1_MB2)         
        selectionMasks['dt_cls_negMB1_MB34'] = buildMask(dt_clusterMasks,DT_sel_negMB1_MB34)         
        #selectionMasks['dt_MCclosure'] = buildMask(dt_clusterMasks,['dr_lep',"dt_jetVeto","dt_muonVeto","dt_time"])         
        selectionMasks['dt_MCclosure'] = buildMask(dt_clusterMasks,['dr_lep',"dt_jetVeto",'dt_muonVeto'])         

        return selectionMasks


    def process(self, events):
        output = self.accumulator.identity()  ## get from histograms
        dataset = events.metadata['dataset']        
                
        isSignal= ak.any(events.gLLP_csc)        
        isData = not(ak.any(events.gLLP_e) or ak.any(events.gLepE))
        if isSignal:        
            csc = ak.any(events.gLLP_csc,axis=1)
            gLLP_dt = ak.firsts(
                ((abs(events.gLLP_decay_vertex_z)<661.0)&\
                 ((events.gLLP_decay_vertex_r<738.0)&(events.gLLP_decay_vertex_r>380.0)))
            )
            gLLP_dt = ak.values_astype(gLLP_dt, np.int32)
            #output['accept'].fill(dataset=dataset,gLLP_csc=ak.firsts(events.gLLP_csc),gLLP_dt=gLLP_dt) ## only 1 LLP
            #events = events[(csc==1)]        

        if isData:
            year = None
            if "_2016" in dataset: year = "2016"
            elif "_2017" in dataset: year = "2017"
            elif "_2018" in dataset: year = "2018" 
            else:
                warnings.warn(" %s does not contain one of the strings: [_2016,_2017,_2018]. No golden json mask applied." % dataset, RuntimeWarning)
            self._year = year 
            if year is not None:           
                events = events[lumiMasks[year](events.runNum,events.lumiSec)]
        output["sumw"][dataset] += len(events)

        ## All possible pairs 
        #cls_lep_pair = ak.cartesian({"cls":cluster_dir,'lep':lep},axis=1,nested=True)
        #dphi_lep_cls = cls_lep_pair.cls.delta_phi(cls_lep_pair.lep)       

        llp      = self.buildLLP(events)
        good_lep,ele,muons = self.buildGoodLeptons(events)
        #reco_muons, glb_muons_pt2, all_muons_pt5 = self.buildRecoMuons(events)
        reco_muons = self.buildRecoMuons(events)
 
        cluster = self.buildCSCcluster(events,good_lep)        

        #if self._saveSkim and not isData:
        #    gParticle = self.buildGenParticles(events)        
        dt_cluster = self.buildDTcluster(events,good_lep)        

        clusterMasks = self.selectCSCcluster(cluster,events) 
        dt_clusterMasks = self.selectDTcluster(dt_cluster,events) 

        #dictionary of cutName:masks
        selectionMasks =   self.buildSelectionMasks(events,good_lep,cluster,clusterMasks,dt_cluster,dt_clusterMasks)
        selectionMasks['loose_ele'] = ak.num( ele[(ele.pt>20) & (abs(ele.eta)<2.5) & (ele.passId)],axis=1)==1
        selectionMasks['loose_mu']  = ak.num( muons[(muons.pt>20) & (abs(muons.eta)<2.4) & (muons.passId)],axis=1)==1

        if self.isElectronChannel:
            preselections = ['trigger_ele','MET',"METfilters",'good_lepton']       
            preselections_loose = ["loose_ele"]       
            if isSignal:
                preselections_csc = preselections + ['cls_PreSel']
                preselections_dt = preselections + ['dt_cls_PreSel']
            else:
                preselections_csc = preselections 
                preselections_dt = preselections 
        else:
            preselections = ['trigger_mu','MET',"METfilters",'good_lepton']       
            preselections_loose = ["loose_mu"]       
            if isSignal:
                preselections_csc = preselections + ['cls_PreSel']
                preselections_dt = preselections + ['dt_cls_PreSel']
            else:
                preselections_csc = preselections 
                preselections_dt = preselections 


        CSC_sel_ABCD = ["dr_lep","ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto","RE12_veto","IntimeCut","timeSpreadCut","ClusterID"]
        DT_sel_ABCD  = ["dr_lep","dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time","dt_deadzones","dt_noise"]
        #DT_sel_ABCD  = ["dr_lep","dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_etaveto","dt_time","dt_deadzones","dt_noise"]

        regions = {
            "PreSel"       :preselections_csc,            
            "JetMuVeto"    :preselections+["cls_JetMuVeto"],
            "JetMuStaVeto" :preselections+["cls_JetMuStaVeto"],
            "ABCD_negME11"      :preselections+["cls_negME11"],
            "ABCD"         :preselections+["cls_ABCD"],            
            "ABCD_OOT"     :preselections+["cls_OOT"],
            "PreSel_dt"    :preselections_dt,
            "ABCD_dt"      :preselections+["dt_cls_ABCD"],            
            "ABCD_dt_MB2"      :preselections+["dt_cls_ABCD_MB2"],            
            "ABCD_dt_MB34"      :preselections+["dt_cls_ABCD_MB34"],            
            "ABCD_dt_OOT"  :preselections+["dt_cls_OOT"],
            "ABCD_dt_OOT_MB2"  :preselections+["dt_cls_OOT_MB2"],
            "ABCD_dt_OOT_MB34"  :preselections+["dt_cls_OOT_MB34"],
            "ABCD_dt_negMB1"       :preselections+["dt_cls_negMB1"],
            "ABCD_dt_negMB1_MB2"       :preselections+["dt_cls_negMB1_MB2"],
            "ABCD_dt_negMB1_MB34"       :preselections+["dt_cls_negMB1_MB34"],
            "JetMuStaVeto_dt" :preselections+["dt_JetMuStaVeto"],
            "MCclosure_dt" :preselections_loose+["dt_MCclosure"],
            "MCclosure" :preselections_loose+["MCclosure"],
            ##"1cls"         :preselections+["n_cls"],            
            #"StatVeto"     :preselections+["cls_StatVeto"],
            #"ele_W_CR"     :['trigger_ele','MET',"METfilters",'good_electron',"W_CR",],
            #"noselection":[],
        }

        #preselection mask
        p = buildMask(selectionMasks, preselections)

        weights = Weights(len(events),storeIndividual=True)
        if not isData:
            corrections.add_pileup_weight(weights, events.npu,'2018')
            corrections.add_Wpt_kfactor(weights, events.gWPt, dataset)
            if self.isElectronChannel:
                corrections.add_electronSFs(weights, ak.firsts(ele),self._year) 
            else:
                corrections.add_muonSFs(weights, ak.firsts(muons),self._year)       
            
            if isSignal and "rwctau" in dataset:
                ## expect dataset = "HNL_*_pl{ctau_old}_rwctau{ctau_new}"
                ctau_old =  float(dataset.split("_")[-2].replace("pl",""))/10     ##ctau in dataset name is in mm
                ctau_new =  float(dataset.split("_")[-1].replace("rwctau",""))/10 ##new ctau needs to be in cm
                corrections.add_ctau_weight(weights, llp.ctau, ctau_old, ctau_new)
            pass

        if self._debug:
            print(dataset)
            print("Weight statistics: %r" % weights.weightStatistics) 
            print("Weight variations: " , weights.variations) 
            print("Weight nominal: " , np.mean(weights.weight())) 

        if self._saveSkim:
            #cut = ak.any(buildMask(selectionMasks,regions["ABCD"]),axis=1)
            #cut = selectionMasks["Acceptance_csc"]
            #fout["MuonSystem"] = {"cluster":cluster[cut],"gParticle":gParticle[cut],"llp":llp[cut],'lep':good_lep[cut]}
            #fout["MuonSystem"] = {"cluster":cluster[cut],"gParticle":gParticle[cut],'lep':good_lep[cut]}
            #fout["MuonSystem"] = {"cluster":cluster[:20],"gParticle":gParticle[:20]}
            #fout.close()
            if self.isElectronChannel: channel ="ele_"
            else: channel ="muon_"
            filename = dataset + "_skim_"+channel + str(time.time()) + ".root"
            #destination = "root://cmseos.fnal.gov//store/user/kkwok/llp/HNL/skim/skim_Mar13/"
            destination = "root://cmseos.fnal.gov//store/user/kkwok/llp/HNL/skim/skim_Jun22/"

            cls_inTime_CR = (p & ak.any(
                                (selectionMasks['cls_ABCD'] |selectionMasks['cls_OOT']) &
                               (cluster.size>200) &
                               (abs(cluster.dphi_cluster_MET)>=0.7),axis=1)
                            )
            dt_cls_inTime_D = (p & ak.any((selectionMasks['dt_cls_ABCD']) &
                               (dt_cluster.size>=150)  &
                               (abs(dt_cluster.dphi_cluster_lep)>=2.8),axis=1)
                            )
            csc_cls_inTime_D = (p & ak.any((selectionMasks['cls_ABCD']) &
                               (cluster.size>=200)  &
                               (abs(cluster.dphi_cluster_lep)>=2.8),axis=1)
                            )

            cut = (dt_cls_inTime_D)|(csc_cls_inTime_D)

            #cut = ak.num(cluster,axis=1)>=1 ## skim testing cut
            if ak.any(cut):
                print("Found events pass skim cut, writing out")
                with uproot.recreate(filename) as fout:
                    #cluster['passABCD'] =  buildMask(selectionMasks,regions["ABCD"])
                    #dt_cluster['passABCD_dt'] =  buildMask(selectionMasks,regions["ABCD_dt"])
                    fout["MuonSystem"] = uproot_writeable(events[cut], events.fields)    # TODO: find out why we can't write all event fields
                    #fout["MuonSystem"] = {"cluster":cluster[cut],"dt_cluster":dt_cluster[cut],"gParticle":gParticle[cut]}
                    #fout["MuonSystem"] = {"runNum":events.runNum[cut],"lumiSec":events.lumiSec[cut],"evtNum":events.evtNum[cut],
                    #                      "nDtRings":events.nDtRings,"
                    #                        "cluster":cluster[cut],"dt_cluster":dt_cluster[cut]}
            
                copyproc = XRootD.client.CopyProcess()
                copyproc.add_job(source = os.path.abspath(os.path.join(".", filename)),
                                 target = (destination + f"/{filename}"), force=True)
                copyproc.prepare()
                print("Copying skim output to  ...",destination )
                copyproc.run()
                client = XRootD.client.FileSystem("root://cmseos.fnal.gov/")
                status = client.locate(
                    "/store/user/kkwok/llp/HNL/skim/"+filename,
                    XRootD.client.flags.OpenFlags.READ,
                )
                assert status[0].ok
                del client
                del copyproc
                
                try:
                    os.remove(os.path.abspath(os.path.join(".", filename)))
                    print("% s removed successfully" % os.path.abspath(os.path.join(".", filename)))
                except OSError as error:
                    print(error)
                    print("File path can not be removed")
            else:
                print("No events pass skim cut")
            return output

        if self._forLimit:
            cutnames = regions["ABCD"]
            cut = buildMask(selectionMasks,cutnames)
            region = "ABCD"
            w_cls      = (weights.weight() * ak.ones_like(cluster.size))[cut] ## use size to pick-up the cluster shape
            output["dphi_cluster_csc"].fill(dataset=dataset,region=region,
                                            ClusterSize=ak.flatten(cluster[cut].size),
                                            dphi_lep =np.abs(ak.flatten(cluster[cut].dphi_cluster_lep)),
                                            dphi_MET=np.abs(ak.flatten(cluster[cut].dphi_cluster_MET)),
                                            weight=ak.flatten(w_cls))
            for region in ["ABCD_dt","ABCD_dt_MB2","ABCD_dt_MB34"]:
                cutnames = regions[region]
                cut = buildMask(selectionMasks,cutnames)
                w_cls      = (weights.weight() * ak.ones_like(dt_cluster.size))[cut] ## use size to pick-up the cluster shape
                output["dphi_cluster_dt"].fill(dataset=dataset,region=region,
                                                ClusterSize=ak.flatten(dt_cluster[cut].size),
                                                dphi_lep =np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_lep)),
                                                dphi_MET=np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_MET)),
                                                weight=ak.flatten(w_cls))
            return output

        if not isData and self._runSys:
            ## CSC systematics

            def fillsys(h,sys_name,cls,cut,weights,isPartial=False):
                ## arrays to be filled
                nhit = ak.flatten(cls[cut].size)
                dphi_lep =np.abs(ak.flatten(cls[cut].dphi_cluster_lep))
                if sys_name=="nominal":
                    ## fill nominal
                    w_cls = ak.flatten((weights.weight() * ak.ones_like(cls.size))[cut]) 
                    h.fill(dataset=dataset,syst=sys_name,ClusterSize=nhit,dphi_lep=dphi_lep,weight=w_cls)
                elif isPartial:
                    ## Partial weight (i.e. all sys but exclude sys_name)
                    w_cls = ak.flatten((weights.partial_weight(exclude=[sys_name]) * ak.ones_like(cls.size))[cut])
                    h.fill(dataset=dataset,syst="no"+sys_name,ClusterSize=nhit,dphi_lep=dphi_lep,weight=w_cls)
                else:
                    ## Varied weight (i.e. all sys but varies sys_name)
                    w_cls = ak.flatten((weights.weight(sys_name) * ak.ones_like(cls.size))[cut])
                    if self._debug:
                        print(sys_name, "weight = ",np.mean(w_cls))
                    h.fill(dataset=dataset,syst=sys_name,ClusterSize=nhit,dphi_lep=dphi_lep,weight=w_cls)

            year = self._year
            if self.isElectronChannel:
                systematics = ["Wpt","WptFSR","WptISR","electron_ID_SF_value","electron_trigger_SF_value"]
            else:
                systematics = ["Wpt","WptFSR","WptISR",f"muon_ID_{year}_value",f"muon_ISO_{year}_value",f"muon_trigger_{year}_value"]
                        
            cut      = buildMask(selectionMasks,regions["ABCD"])
            cut_dt   = buildMask(selectionMasks,regions["ABCD_dt"])
            cut_dt_mb2   = buildMask(selectionMasks,regions["ABCD_dt_MB2"])
            cut_dt_mb34   = buildMask(selectionMasks,regions["ABCD_dt_MB34"])
            ## Nominal 
            fillsys(output['dphi_cluster_syst']   ,"nominal",    cluster,cut   ,weights)
            fillsys(output['dphi_cluster_dt_syst'],"nominal"    , dt_cluster,cut_dt,weights)
            fillsys(output['dphi_cluster_dt_mb2_syst'],"nominal", dt_cluster,cut_dt_mb2,weights)
            fillsys(output['dphi_cluster_dt_mb34_syst'],"nominal",dt_cluster,cut_dt_mb34,weights)
            ## Exclude Wpt to see the effect of WpT reweighting
            fillsys(output['dphi_cluster_syst']   ,"Wpt",    cluster,cut   ,weights,isPartial=True)
            fillsys(output['dphi_cluster_dt_syst'],"Wpt", dt_cluster,cut_dt,weights,isPartial=True)
            ## Vaiations 
            for sys in systematics:
                for suffix in ["Up","Down"]:
                    sys_name = sys+suffix
                    fillsys(output['dphi_cluster_syst']   ,sys_name,    cluster,cut   ,weights)
                    fillsys(output['dphi_cluster_dt_syst'],sys_name, dt_cluster,cut_dt,weights)
            ## llp eta:
            cut = ak.any(cut,axis=1)
            cut_dt = ak.any(cut_dt,axis=1)

            output['gLLP_eta_csc'].fill(dataset=dataset,gllp_eta = abs(ak.firsts(llp[cut].eta)), weight=weights.weight()[cut])
            output['gLLP_eta_dt'].fill(dataset=dataset ,gllp_eta = abs(ak.firsts(llp[cut_dt].eta)), weight=weights.weight()[cut_dt])
            output['gLLP_phi'].fill(dataset=dataset,region= "ABCD"   , gLLP_phi = abs(ak.firsts(llp[cut].phi)), weight=weights.weight()[cut])
            output['gLLP_phi'].fill(dataset=dataset,region= "ABCD_dt",gLLP_phi = abs(ak.firsts(llp[cut_dt].phi)), weight=weights.weight()[cut_dt])
            return output

            
        ## Fill no selection plots
        output['nLeptons'].fill(dataset=dataset, nLeptons = events.nLeptons, weight=weights.weight())
        
        output['elePt'].fill(dataset=dataset , elePt  = ak.flatten(ele.pt) )
        output['eleEta'].fill(dataset=dataset, eleEta = ak.flatten(ele.eta))
        output['muPt'].fill(dataset=dataset , muPt  = ak.flatten(muons.pt) )
        output['muEta'].fill(dataset=dataset, muEta = ak.flatten(muons.eta))

        output['nJets'].fill(dataset=dataset, nJets = events.nJets, weight=weights.weight())
        output['jetMet_dPhi'].fill(dataset=dataset, jetMet_dPhi = events.jetMet_dPhi, weight=weights.weight())

        output["nPU"].fill(dataset=dataset,nPU=events.npu,weight=weights.weight())        
        output["gWPt"].fill(dataset=dataset,gWPt=events.gWPt,weight=weights.weight())        
        output["gWPt_noweight"].fill(dataset=dataset,gWPt=events.gWPt)        
        output["nPU_noweight"].fill(dataset=dataset,nPU=events.npu)        
        if isSignal:
            output['accept'].fill(dataset=dataset,
                                  gLLP_csc=ak.firsts(events.gLLP_csc),
                                  gLLP_dt=gLLP_dt,weight=weights.weight()) ## only 1 LLP
            cut = selectionMasks["Acceptance_csc"]
            output['gLLP_e'].fill(dataset=dataset  ,region="gLLP_csc" ,gLLP_e = ak.firsts(llp[cut].e) , weight=weights.weight()[cut])
            output['gLLP_pt'].fill(dataset=dataset ,region="gLLP_csc" ,gLLP_pt = ak.firsts(llp[cut].pt), weight=weights.weight()[cut])
            output['gLLP_eta'].fill(dataset=dataset,region="gLLP_csc" ,gLLP_eta = ak.firsts(llp[cut].eta), weight=weights.weight()[cut])
            output['glepdPhi'].fill(dataset=dataset,region="gLLP_csc" ,gLLP_lepdPhi = np.abs(ak.flatten(events[cut].gLLP_lepdPhi)), weight=weights.weight()[cut])
            output["metXYCorr"].fill(dataset=dataset,region="gLLP_csc",metXYCorr=events[cut].metXYCorr,weight=weights.weight()[cut]) 
            cut = selectionMasks["Acceptance_dt"]
            output['gLLP_e'].fill(dataset=dataset  ,region="gLLP_dt" ,gLLP_e = ak.firsts(llp[cut].e) , weight=weights.weight()[cut])
            output['gLLP_pt'].fill(dataset=dataset ,region="gLLP_dt" ,gLLP_pt = ak.firsts(llp[cut].pt), weight=weights.weight()[cut])
            output['gLLP_eta'].fill(dataset=dataset,region="gLLP_dt" ,gLLP_eta = ak.firsts(llp[cut].eta), weight=weights.weight()[cut])
            output['glepdPhi'].fill(dataset=dataset,region="gLLP_dt" ,gLLP_lepdPhi = np.abs(ak.flatten(events[cut].gLLP_lepdPhi)), weight=weights.weight()[cut])
            output["metXYCorr"].fill(dataset=dataset,region="gLLP_dt",metXYCorr=events[cut].metXYCorr,weight=weights.weight()[cut]) 

            ## get CSC cluster masks

            #Events with clusterID pass
            #llp_selection = maskAndFill(llp.e,ak.any(cluster[cut].llp_match,axis=1),len(llp.e[0])*[0])
            #Events with any cluster matching to llp
            #llp_selection = ak.values_astype( ak.any(cluster.llp_match,axis=1),np.int32 )
            llp_dau = util.pack(events,"gLLP_daughter_")
            for region in [50,100,130,150,200]:
                br = ak.any(abs(llp_dau.id)<=4,axis=1)
                had_e = ak.sum(ak.mask(llp_dau.e,(abs(llp_dau.id)<=4)),axis=1)
                cut = (selectionMasks["Acceptance_csc"])&(br)
                llp_selection = ak.values_astype(( (ak.num(cluster,axis=1)==1)&(br) & (ak.all((cluster.llp_match) & (cluster.size>=region), axis=1))),np.int32)
                #llp_selection = ak.values_astype((ak.num(cluster,axis=1)==1 & (ak.all((cluster.llp_match) & (cluster.size>=region), axis=1))),np.int32)
                #llp_selection = ak.values_astype( ak.any(cluster.llp_match & (cluster.dr_cluster_lep>0.8),axis=1),np.int32 )

                output['llp_cls_eff_z'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],z=ak.flatten(abs(llp.z[cut])),weight=weights.weight()[cut])
                output['llp_cls_eff_r'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],r=ak.flatten(llp.r[cut]),weight=weights.weight()[cut])
                output['llp_cls_eff_e'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],e=ak.flatten(llp.e[cut]),weight=weights.weight()[cut])
                #output['llp_cls_eff_e'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],e=had_e,weight=weights.weight()[cut])

                cut = (selectionMasks["Acceptance_dt"])&(br)
                #llp_selection = ak.values_astype( ak.any(dt_cluster.llp_match,axis=1),np.int32 )
                llp_selection = ak.values_astype(((ak.num(dt_cluster,axis=1)==1)&(br) & (ak.all((dt_cluster.llp_match) & (dt_cluster.size>=region), axis=1))),np.int32)
                #llp_selection = ak.values_astype((ak.num(dt_cluster,axis=1)==1 & (ak.all((dt_cluster.llp_match) & (dt_cluster.size>=region), axis=1))),np.int32)
                #llp_selection = ak.values_astype( ak.any(dt_cluster.llp_match & (dt_cluster.dr_cluster_lep>0.8),axis=1),np.int32 )
                output['llp_cls_dt_eff_z'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],z=ak.flatten(abs(llp.z[cut])),weight=weights.weight()[cut])
                output['llp_cls_dt_eff_r'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],r=ak.flatten(llp.r[cut]),weight=weights.weight()[cut])
                output['llp_cls_dt_eff_e'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],e=ak.flatten(llp.e[cut]),weight=weights.weight()[cut])
                #output['llp_cls_dt_eff_e'].fill(dataset=dataset,region=str(region),selection=llp_selection[cut],e=had_e,weight=weights.weight()[cut])


        output["metXYCorr"].fill(dataset=dataset,region="noselection",metXYCorr=events.metXYCorr,weight=weights.weight()) 

        for region in ["CleanLep","CleanLep_Tight"]:
            if region=="CleanLep":
                acc_csc = "Acceptance_csc"
                acc_dt = "Acceptance_dt"
            else:
                acc_csc = "Acceptance_csc_tight"
                acc_dt = "Acceptance_dt_tight"
            if isSignal:
                #preselection_mask =  buildMask(selectionMasks,preselections+["Acceptance_csc"])
                preselection_mask =  buildMask(selectionMasks,preselections+[acc_csc])
                cut = (preselection_mask) & (cluster.dr_cluster_lep>0.8)&(cluster.time < 12.5)&(cluster.time>-5)&(cluster.llp_match)
            else:
                preselection_mask =  buildMask(selectionMasks,preselections)
                cut = (preselection_mask) & (cluster.dr_cluster_lep>0.8)&(cluster.time < 12.5)&(cluster.time>-5)
    
            w_cls      = (weights.weight() * ak.ones_like(cluster.size))[cut] ## use size to pick-up the cluster shape
            output["ClusterJetVetoPt"].fill(dataset=dataset,region=region,ClusterJetVetoPt=ak.flatten(cluster[cut].JetVetoPt),weight=ak.flatten(w_cls))        
            output["ClusterMuonVetoPt"].fill(dataset=dataset,region=region,ClusterMuonVetoPt=ak.flatten(cluster[cut].MuonVetoPt),weight=ak.flatten(w_cls))              
            output["ClusterME11_12"].fill(dataset=dataset,region=region,ClusterME11_12=ak.flatten(cluster[cut].ME11_12),weight=ak.flatten(w_cls))        
            if isSignal:
                preselection_mask =  buildMask(selectionMasks,preselections+[acc_dt])
                cut = (preselection_mask) & (dt_cluster.dr_cluster_lep>0.8)&(dt_cluster.rpcBx==0)&(dt_cluster.llp_match)
            else:
                preselection_mask =  buildMask(selectionMasks,preselections)
                #cut = (preselection_mask) & (dt_cluster.dr_cluster_lep>0.8)&(dt_cluster.rpcBx==0)&(dt_cluster.MuonVetoLooseId)
                cut = (preselection_mask) & (dt_cluster.dr_cluster_lep>0.8)&(dt_cluster.rpcBx==0)
            w_cls      = (weights.weight() * ak.ones_like(dt_cluster.size))[cut] ## use size to pick-up the cluster shape
            output["ClusterJetVetoPt_dt"].fill(dataset=dataset,region=region,ClusterJetVetoPt=ak.flatten(dt_cluster[cut].JetVetoPt),weight=ak.flatten(w_cls))        
            output["ClusterMuonVetoPt_dt"].fill(dataset=dataset,region=region,ClusterMuonVetoPt=ak.flatten(dt_cluster[cut].MuonVetoPt),weight=ak.flatten(w_cls))       
            output["ClusterMuonVetoPt_dt"].fill(dataset=dataset,region=region,ClusterMuonVetoPt=ak.flatten(dt_cluster[cut].MuonVetoPt),weight=ak.flatten(w_cls))       
            output["ClusterMB1_dt"].fill(dataset=dataset,region=region,ClusterMB1=ak.flatten(dt_cluster[cut].nMB1),weight=ak.flatten(w_cls))        

 
        ## Fill regions plot
        for region,cutnames in regions.items():

            ## Fill other regions without cutflows
            cut = buildMask(selectionMasks,cutnames)
            if self._debug:
                if "MC" in region: print(region," cuts = ", cutnames)
            if cut.ndim==1:
                ev_cut = cut                  ##This is a per-event cut
            else:
                ev_cut = ak.any(cut,axis=1)   ##This is a per-cluster cut, require at least 1 passing cluster
 
            ev_cut = ak.fill_none(ev_cut,False)
            w_evt = weights.weight()[ev_cut]

            if not "dt" in region:
                w_cls      = (weights.weight() * ak.ones_like(cluster.size))[cut] ## use size to pick-up the cluster shape
                output["dphi_cluster_csc"].fill(dataset=dataset,region=region,
                                                ClusterSize=ak.flatten(cluster[cut].size),
                                                dphi_lep =np.abs(ak.flatten(cluster[cut].dphi_cluster_lep)),
                                                dphi_MET=np.abs(ak.flatten(cluster[cut].dphi_cluster_MET)),
                                                weight=ak.flatten(w_cls))
                output["ClusterSize"].fill(dataset=dataset,region=region,
                                           ClusterSize=ak.flatten(cluster[cut].size),
                                           weight=ak.flatten(w_cls))        
                output["ClusterTime"].fill(dataset=dataset,region=region,
                                           ClusterTime=ak.flatten(cluster[cut].time),
                                           weight=ak.flatten(w_cls))        
                output["ClusterTimeSpread"].fill(dataset=dataset,region=region,
                                           ClusterTimeSpread=ak.flatten(cluster[cut].timeSpread),
                                           weight=ak.flatten(w_cls))        
                output["ClusterEta"].fill(dataset=dataset,region=region,
                                           ClusterEta=np.abs(ak.flatten(cluster[cut].eta)),
                                           weight=ak.flatten(w_cls))        
                output["ClusterAvgStation10"].fill(dataset=dataset,region=region,
                                           ClusterAvgStation10=np.abs(ak.flatten(cluster[cut].AvgStation10)),
                                           weight=ak.flatten(w_cls))        
                output["ClusterNStation10"].fill(dataset=dataset,region=region,
                                           ClusterNStation10=ak.flatten(cluster[cut].NStation10),
                                           weight=ak.flatten(w_cls))        
                output["ClusterME11_12"].fill(dataset=dataset,region=region,
                                           ClusterME11_12=ak.flatten(cluster[cut].ME11_12),
                                           weight=ak.flatten(w_cls))        
                output["ClusterJetVetoPt"].fill(dataset=dataset,region=region,
                                           ClusterJetVetoPt=ak.flatten(cluster[cut].JetVetoPt),
                                           weight=ak.flatten(w_cls))        
                output["ClusterMuonVetoPt"].fill(dataset=dataset,region=region,
                                           ClusterMuonVetoPt=ak.flatten(cluster[cut].MuonVetoPt),
                                           weight=ak.flatten(w_cls))        
                output["metXYCorr"].fill(dataset=dataset,region=region,
                                         metXYCorr=events[ev_cut].metXYCorr,
                                        weight=w_evt)                    
                output['jetPt'].fill(dataset=dataset, region = region, 
                                    jetPt = ak.to_numpy(ak.firsts(events[ev_cut].jetPt)),
                                    weight=w_evt)
                output["MT"].fill(dataset=dataset,region=region,MT=events[ev_cut].MT,weight=w_evt)       
            else:
                w_cls      = (weights.weight() * ak.ones_like(dt_cluster.size))[cut] ## use size to pick-up the cluster shape

                if region=="ABCD_dt_OOT" or region=="ABCD_dt":
                    if isData:
                        output["Cluster_runNum_dt"].fill(dataset=dataset,region=region,ClusterSize=ak.flatten(dt_cluster[cut].size),
                                                RunNumber = ak.flatten((ak.ones_like(dt_cluster.size)[cut])*(events.runNum)),
                                                weight=ak.flatten(w_cls))
                        output["Cluster_rz_dt"].fill(dataset=dataset,region=region,R=ak.flatten(dt_cluster[cut].R),
                                                Z = ak.flatten(dt_cluster[cut].Z),
                                                weight=ak.flatten(w_cls))
                        output["Cluster_phi_dt"].fill(dataset=dataset,region=region,phi=ak.flatten(dt_cluster[cut].phi),
                                                weight=ak.flatten(w_cls))


                output["dphi_cluster_dt"].fill(dataset=dataset,region=region,
                                                ClusterSize=ak.flatten(dt_cluster[cut].size),
                                                dphi_lep =np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_lep)),
                                                dphi_MET=np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_MET)),
                                                weight=ak.flatten(w_cls))
                output["ClusterSize_dt"].fill(dataset=dataset,region=region,
                                           ClusterSize=ak.flatten(dt_cluster[cut].size),
                                           weight=ak.flatten(w_cls))        
                output["ClusterTime_dt"].fill(dataset=dataset,region=region,
                                           ClusterBx=ak.flatten(dt_cluster[cut].rpcBx),
                                           weight=ak.flatten(w_cls))        
                output["ClusterEta_dt"].fill(dataset=dataset,region=region,
                                           ClusterEta=np.abs(ak.flatten(dt_cluster[cut].eta)),
                                           weight=ak.flatten(w_cls))        
                output["ClusterAvgStation10_dt"].fill(dataset=dataset,region=region,
                                           ClusterAvgStation10=np.abs(ak.flatten(dt_cluster[cut].AvgStation10)),
                                           weight=ak.flatten(w_cls))        
                output["ClusterNStation10_dt"].fill(dataset=dataset,region=region,
                                           ClusterNStation10=ak.flatten(dt_cluster[cut].NStation10),
                                           weight=ak.flatten(w_cls))       
                output["ClusterMB1_dt"].fill(dataset=dataset,region=region,
                                           ClusterMB1=ak.flatten(dt_cluster[cut].nMB1),
                                           weight=ak.flatten(w_cls))        
                output["ClusterJetVetoPt_dt"].fill(dataset=dataset,region=region,
                                           ClusterJetVetoPt=ak.flatten(dt_cluster[cut].JetVetoPt),
                                           weight=ak.flatten(w_cls))        
                output["ClusterMuonVetoPt_dt"].fill(dataset=dataset,region=region,
                                           ClusterMuonVetoPt=ak.flatten(dt_cluster[cut].MuonVetoPt),
                                           weight=ak.flatten(w_cls))       
           
 
        ## fill cutflow plots:
        output['cutflow'].fill(dataset=dataset,region="csc_cutflow",cutflow="NoSelection",weight=weights.weight())
        output['cutflow'].fill(dataset=dataset,region="dt_cutflow",cutflow="NoSelection",weight=weights.weight())
        if isSignal:
            acc_csc_preselections = ["Acceptance_csc"] + preselections + ["n_cls"]
            acc_dt_preselections  = ["Acceptance_dt" ] + preselections + ["n_cls_dt"]
        else:
            acc_csc_preselections =  preselections+ ["n_cls"]   
            acc_dt_preselections  =  preselections+ ["n_cls_dt"]

        for i,sel in enumerate(acc_csc_preselections):
            allcuts= buildMask(selectionMasks,acc_csc_preselections[0:i+1])       ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="csc_cutflow",cutflow=sel,weight=weights.weight()[allcuts])
        for i,sel in enumerate(acc_dt_preselections):
            allcuts= buildMask(selectionMasks,acc_dt_preselections[0:i+1])       ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow",cutflow=sel,weight=weights.weight()[allcuts])

        allPreSel= buildMask(selectionMasks,acc_csc_preselections) 
        allPreSel_dt= buildMask(selectionMasks,acc_dt_preselections) 
        
        for i,sel in enumerate(CSC_sel_ABCD):
            allcuts= ak.any( (buildMask(clusterMasks,CSC_sel_ABCD[0:i+1]) & allPreSel), axis=1)     ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="csc_cutflow",cutflow=sel,weight=weights.weight()[allcuts])
        for i,sel in enumerate(DT_sel_ABCD):
            allcuts= ak.any( (buildMask(dt_clusterMasks,DT_sel_ABCD[0:i+1]) & allPreSel_dt), axis=1)     ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow",cutflow=sel,weight=weights.weight()[allcuts])


        return output

    def postprocess(self, accumulator):
        # set everything to 1/fb scale
        lumi = 1000  # [1/pb]

        scale = {}
        for dataset, dataset_sumw in accumulator['sumw'].items():
            if "rwctau" in dataset:
                ctau    = float(dataset.split("_")[-1].replace("rwctau",""))
                mass    = float(dataset.split("_")[2].replace("mHNL","").replace("p","."))
                if "tau" in dataset:
                    src_sample       = "_".join(dataset.split("_")[:-1])
                    src_ctau         = float(src_sample.split("_")[-1].replace("pl","")) 
                    src_xsecTimeseff = corrections.load_xsection()[src_sample]  ## with eff
                    src_xsec         = corrections.reweightXsec(src_ctau,mass,True)
                    filter_eff       =  src_xsecTimeseff/src_xsec
                    xsec             = corrections.reweightXsec(ctau,mass,True)*filter_eff
                else:
                    xsec    = corrections.reweightXsec(ctau,mass)
                if self._debug:
                    print("reweighting ct for ",dataset," with xsec = ",xsec)
                scale[dataset] = lumi*xsec/dataset_sumw
                
            elif dataset in corrections.load_xsection().keys():
                scale[dataset] = lumi*corrections.load_xsection()[dataset]/dataset_sumw
                if self._debug:
                    print(" dataset =",dataset,"  xsec = ",corrections.load_xsection()[dataset])
            else:
                warnings.warn("Missing cross section for dataset %s.  No normalization applied. " % dataset, RuntimeWarning)
                #scale[dataset] = lumi / dataset_sumw

        if self._debug:
            print("Scaling with scale = " , scale)
        for h in accumulator.values():
            if isinstance(h, hist.Hist):
                h.scale(scale, axis="dataset")

        return accumulator

