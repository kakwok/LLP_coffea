import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
from coffea import hist, processor
import uproot
# register our candidate behaviors
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema, TreeMakerSchema
import time

from coffea.analysis_tools import Weights, PackedSelection

import HNLprocessor.corrections as corrections 
from  HNLprocessor.histograms import histograms 
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
    def __init__(self,isElectronChannel=True,saveSkim=False,debug=False):
        self._debug = debug
        self._saveSkim = saveSkim
        self.isElectronChannel = isElectronChannel
        self.isMuonChannel = not(isElectronChannel)
        ##define histograms 
        histograms['sumw']= processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator( histograms )

    @property
    def accumulator(self):
        return self._accumulator

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
        good_ele = ele[(ele.pt>35) & (abs(ele.eta)<2.4) & (ele.passId)]
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
            'z':events.gLLP_decay_vertex_z ,
            'r':events.gLLP_decay_vertex_r,
            'ctau':events.gLLP_ctau,
        }) 
        return llp

    def buildRecoMuon(self,events):
        muons= ak.zip({
            'pt': events.muonPt,
            'e': events.muonE,
            'eta': events.muonEta,
            'phi': events.muonPhi,
            'Iso': events.muonIso,
            'IsGlobal': events.muonIsGlobal,
            'TightId': events.muonTightId,
            'LooseId': events.muonLooseId,
        })
        return muons

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
                "Cluster_match_gParticle_id":events.cscRechitCluster3_match_gParticle_id,
                "Cluster_match_gParticle_minDeltaR":events.cscRechitCluster3_match_gParticle_minDeltaR,
                "Cluster_match_gParticle":events.cscRechitCluster3_match_gParticle,
                "Cluster_match_gParticle_index":events.cscRechitCluster3_match_gParticle_index,
                "Cluster_match_gParticle_eta":events.cscRechitCluster3_match_gParticle_eta,
                "Cluster_match_gParticle_phi":events.cscRechitCluster3_match_gParticle_phi,
                "Cluster_match_gParticle_e":events.cscRechitCluster3_match_gParticle_E,
                "Cluster_match_gParticle_pt":events.cscRechitCluster3_match_gParticle_pt,
            }
        )
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
        
        muonVeto_old = ~((muVeto.pt>20) & abs(muVeto.eta<2.4)) 
        if self.isMuonChannel: muonVeto = muonVeto_old & (cluster.dr_cluster_lep>0.8)
        jetVeto  = ~((jetVeto.pt>10)& abs(jetVeto.eta<2.4))
        RE12_veto     = (cluster.RE12==0)
        MB1seg_veto   = (cluster.MB1seg==0)
        RB1_veto      = (cluster.RB1==0)
        ME11_12_veto  = (cluster.ME11_12==0)
        ME11_12_veto_inv = ~(cluster.ME11_12==0)

        oot_timecut   = (cluster.time < -12.5) #OOT for data
        IntimeCut     = (cluster.time < 12.5) & (cluster.time>-5) ## In-time otherwise        
        timeSpreadCut = (cluster.timeSpread<20)        
        dphi_met      = (abs(cluster.dphi_cluster_MET)<0.75)        
        dphi_lep      = (abs(cluster.dphi_cluster_lep)>2.5)      
        dr_lep      = (cluster.dr_cluster_lep>0.4)
        
        pid_gamma     = (cluster.Cluster_match_gParticle_id==22)
        pid_gluon     = (cluster.Cluster_match_gParticle_id==21)
        pid_quarks    = (abs(cluster.Cluster_match_gParticle_id)<10)
        pid_muons     = (abs(cluster.Cluster_match_gParticle_id)==13) & (cluster.Cluster_match_gParticle_minDeltaR<0.4) 
        pid_mesons    = (abs(cluster.Cluster_match_gParticle_id)>100) & (abs(cluster.Cluster_match_gParticle_id)<1000)
        
        NRechit100 = (cluster.size<100)
        
        clusterMasks = ak.zip({
            "ClusterID"     : ClusterID     ,  
            "muonVeto" : muonVeto ,
            "jetVeto"  : jetVeto  ,
            "RE12_veto"     : RE12_veto     ,
            "MB1seg_veto"   : MB1seg_veto   ,
            "RB1_veto"      : RB1_veto      ,
            "ME11_12_veto"  : ME11_12_veto  ,
            "ME11_12_veto_inv"  : ME11_12_veto_inv  ,                             
            "OOT_timeCut"   : oot_timecut   ,
            "IntimeCut"     : IntimeCut     ,
            "timeSpreadCut" : timeSpreadCut ,
            "dphi_MET"      : dphi_met      ,
            "dphi_lep"      : dphi_lep      ,
            "dr_lep"        : dr_lep        ,
            "pid_gamma"     : pid_gamma     ,
            "pid_gluon"     : pid_gluon     ,
            "pid_quarks"    : pid_quarks    ,
            "genMuon"       : pid_muons     ,
            "pid_mesons"    : pid_mesons    ,
            "muonVeto_old"  : muonVeto_old,
            "NRechit100"    : NRechit100    ,
            })
        return clusterMasks 


    def buildDTcluster(self,events,good_lep):
        
        dt_cluster = ak.zip(
            {
                 "size":events.dtRechitClusterSize,
                 "x":events.dtRechitClusterX,
                 "y":events.dtRechitClusterY,
                 "z":events.dtRechitClusterZ,
                 "llp_x":events.dtRechitCluster_match_gLLP_decay_x,
                 "llp_y":events.dtRechitCluster_match_gLLP_decay_y,
                 "llp_z":events.dtRechitCluster_match_gLLP_decay_z,
                 "llp_match":events.dtRechitCluster_match_gLLP,
                 "eta":events.dtRechitClusterEta,
                 "phi":events.dtRechitClusterPhi,
                 "JetVetoPt":events.dtRechitClusterJetVetoPt,
                 "MuonVetoPt":events.dtRechitClusterMuonVetoPt,
                 "MuonVetoEta":events.dtRechitClusterMuonVetoEta,
                 "MuonVetoPhi":events.dtRechitClusterMuonVetoPhi,
                 "MuonVetoE":events.dtRechitClusterMuonVetoE,
                 #"MuonVetoTightId":events.dtRechitClusterMuonVetoTightId,
                 #"MuonVetoTightIso":events.dtRechitClusterMuonVetoTightIso,
                 #"MuonVetoLooseId":events.dtRechitClusterMuonVetoLooseId,
                 "NStation10":events.dtRechitClusterNStation10,
                 "AvgStation10":events.dtRechitClusterAvgStation10,
                 "MaxStation":events.dtRechitClusterMaxStation,
                 "nRPC":events.dtRechitCluster_match_RPChits_dPhi0p5,
                 "nMB1_cosmic_minus":events.dtRechitCluster_match_MB1hits_cosmics_minus,
                 "nMB1_cosmic_plus":events.dtRechitCluster_match_MB1hits_cosmics_plus,
                 "nMB1":events.dtRechitCluster_match_MB1hits_0p5,
                 "rpcBx":events.dtRechitCluster_match_RPCBx_dPhi0p5,
                 "dphi_cluster_MET":events.dtRechitClusterMetEENoise_dPhi, 
                 "Cluster_match_gParticle_id":events.dtRechitCluster_match_gParticle_Id,
                 "Cluster_match_gParticle_minDeltaR":events.dtRechitCluster_match_gParticle_deltaR, 
                 "Cluster_match_gParticle_eta":events.dtRechitCluster_match_gParticle_Eta,
                 "Cluster_match_gParticle_phi":events.dtRechitCluster_match_gParticle_Phi,
                 "Cluster_match_gParticle_e":events.dtRechitCluster_match_gParticle_E,
                 "Cluster_match_gParticle_pt":events.dtRechitCluster_match_gParticle_Pt,
                 "NSegmentStation1":events.dtRechitClusterNSegmentStation1,
                 "NSegmentStation2":events.dtRechitClusterNSegmentStation2,
                 "NSegmentStation3":events.dtRechitClusterNSegmentStation3,
                 "NSegmentStation4":events.dtRechitClusterNSegmentStation4,
                 "MaxStation" : events.dtRechitClusterMaxStation,
                 "MaxStationRatio" : events.dtRechitClusterMaxStationRatio,
                 "NChamber":events.dtRechitClusterNChamber,
                 "MaxChamber":events.dtRechitClusterMaxChamber,
                 },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )

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
        
        dphi_dt_cluster_lep = ak.fill_none(dt_cluster.delta_phi(ak.firsts(good_lep)),-999,axis=None)
        dr_dt_cluster_lep = ak.fill_none(dt_cluster.delta_r(ak.firsts(good_lep)),-999,axis=None)
        
        dt_cluster = ak.with_field(dt_cluster,dphi_dt_cluster_lep,"dphi_cluster_lep")
        dt_cluster = ak.with_field(dt_cluster,dr_dt_cluster_lep,"dr_cluster_lep")
        
        """
        muons= ak.zip(
            {
                'pt': events.muonPt,
                'e': events.muonE,
                'eta': events.muonEta,
                'phi': events.muonPhi,
                'Iso': events.muonIso,
                'IsGlobal': events.muonIsGlobal,
                'TightId': events.muonTightId,
                'LooseId': events.muonLooseId,
            }
        )
        cls_muons_pairs = ak.cartesian({"cls":dt_cluster_dir,
                                        "muon":muons[muons.pt>10]},axis=1,nested=True)
        delta_r_muon = cls_muons_pairs.cls.delta_r(cls_muons_pairs.muon)


        dt_cluster = ak.with_field(dt_cluster,ak.any(delta_r_muon<0.4,axis=2),"MuonVeto_0p4")
        dt_cluster = ak.with_field(dt_cluster,ak.any(delta_r_muon<0.8,axis=2),"MuonVeto_0p8")
        """
        return dt_cluster
 
    def selectDTcluster(self,dt_cluster,events):
        dt_jetVeto  = (dt_cluster.JetVetoPt<20.0)
        dt_jetVeto_inv = (dt_cluster.JetVetoPt>20.0)
        dt_muonVeto = (dt_cluster.MuonVetoPt<10.0)
        
        
        dt_MB1veto  = (dt_cluster.nMB1<=1)
        dt_MB1veto_inv = (dt_cluster.nMB1>1)
        
        dt_RPC      = (dt_cluster.nRPC>=1)
        dt_MB1adj   = (dt_cluster.nMB1_cosmic_minus<=8) & (dt_cluster.nMB1_cosmic_plus<=8)
        dt_time     = (dt_cluster.rpcBx==0)
        dt_OOT      = (dt_cluster.rpcBx>=-100)&(dt_cluster.rpcBx<0)
        dt_dphi_MET  = (abs(dt_cluster.dphi_cluster_MET)<0.75)
        dt_size      = (dt_cluster.size>=100)
        dr_lep      = (dt_cluster.dr_cluster_lep>0.4)
        #MuonVeto_0p4 = ~(dt_cluster.MuonVeto_0p4)
        #MuonVeto_0p8 = ~(dt_cluster.MuonVeto_0p8)
        dt_rechits_130 = (dt_cluster.size>=130)
        dt_dphi_lep_3p0 = (dt_cluster.dphi_cluster_lep>=3.0)
        
        pid_muons     = (abs(dt_cluster.Cluster_match_gParticle_id)==13) & (dt_cluster.Cluster_match_gParticle_minDeltaR<0.4)
        
        dt_deadzones = ~(dt_cluster.Deadzone_1) & ~(dt_cluster.Deadzone_2) 

        segmentsStation_veto = ~(((dt_cluster.MaxStation==2) & (dt_cluster.NSegmentStation1 ==1)) |\
                                ((dt_cluster.MaxStation==3) & ((dt_cluster.NSegmentStation2 ==1) | (dt_cluster.NSegmentStation1 ==1))) |\
                                ((dt_cluster.MaxStation==4) & ((dt_cluster.NSegmentStation3 ==1) | (dt_cluster.NSegmentStation2 ==1) | (dt_cluster.NSegmentStation1 ==1)))) 
        NRechit70 = (dt_cluster.size<70)
        clusterMasks = ak.zip({
                "dt_jetVeto"  :dt_jetVeto  ,
                "dt_muonVeto" :dt_muonVeto ,
                "dt_MB1veto"  :dt_MB1veto  ,
                "dt_RPC"      :dt_RPC      ,
                "dt_MB1adj"   :dt_MB1adj   ,
                "dt_time"     :dt_time     ,
                "dt_OOT"      :dt_OOT      ,
                "dt_dphi_MET" :dt_dphi_MET ,
                "dt_size"     :dt_size     ,
                "dr_lep"      :dr_lep      ,
                "dt_jetVeto_inv" : dt_jetVeto_inv ,
                "dt_MB1veto_inv" : dt_MB1veto_inv ,
                "dt_rechits_130" : dt_rechits_130,
                "dt_dphi_lep_3p0": dt_dphi_lep_3p0,
                "dt_genMuon" : pid_muons,
                "dt_deadzones": dt_deadzones, 
                "segmentsStation_veto" : segmentsStation_veto,
                "NRechit70" : NRechit70,
                #"MuonVeto_0p4":MuonVeto_0p4,
                #"MuonVeto_0p8":MuonVeto_0p8,
                })
        return clusterMasks
    
    def process(self, events):
        output = self.accumulator.identity()  ## get from histograms
        dataset = events.metadata['dataset']        
        
        #start,stop = events._branchargs['entry_start'],events._branchargs['entry_stop']
        #events = uproot.lazy(events._tree)
        #events = events[start:stop]
                
        isSignal= ak.any(events.gLLP_csc)        
        isData = not(ak.any(events.gLLP_e) or ak.any(events.gLepE))
        output["sumw"][dataset] += len(events)
        if isSignal:        
            csc = ak.any(events.gLLP_csc,axis=1)
            gLLP_dt = ak.firsts(
                ((abs(events.gLLP_decay_vertex_z)<661.0)&\
                 ((events.gLLP_decay_vertex_r<738.0)&(events.gLLP_decay_vertex_r>380.0)))
            )
            gLLP_dt = ak.values_astype(gLLP_dt, np.int)
            #output['accept'].fill(dataset=dataset,gLLP_csc=ak.firsts(events.gLLP_csc),gLLP_dt=gLLP_dt) ## only 1 LLP
            #events = events[(csc==1)]        
                                
        
        def buildMask(allMasks,cutnames):
            allcuts = allMasks[cutnames[0]]
            for i,cutname in enumerate(cutnames):
                allcuts = (allcuts) & allMasks[cutname]
            return allcuts





       ## All possible pairs 
        #cls_lep_pair = ak.cartesian({"cls":cluster_dir,'lep':lep},axis=1,nested=True)
        #dphi_lep_cls = cls_lep_pair.cls.delta_phi(cls_lep_pair.lep)       

        llp      = self.buildLLP(events)
        good_lep,ele,muons = self.buildGoodLeptons(events) 
        cluster = self.buildCSCcluster(events,good_lep)        
        if self._saveSkim:
            gParticle = self.buildGenParticles(events)        
        dt_cluster = self.buildDTcluster(events,good_lep)        

        clusterMasks = self.selectCSCcluster(cluster,events) 
        dt_clusterMasks = self.selectDTcluster(dt_cluster,events) 

        #dictionary of cutName:masks
        selectionMasks =   {}

        selectionMasks['Acceptance_csc']   =ak.firsts(events.gLLP_csc)==1
        selectionMasks['Acceptance_dt']=ak.firsts(events.gLLP_dt)==1
        selectionMasks['METfilters']   =events.Flag2_all==True
        selectionMasks['trigger_ele']  =events.SingleEleTrigger==True
        selectionMasks['trigger_mu']   =events.SingleMuonTrigger==True
        selectionMasks['good_lepton']  =ak.num(good_lep,axis=1)==1
        selectionMasks['MET']          =events.metEENoise>=30
        selectionMasks['n_cls']        =ak.num(cluster,axis=1)>=1
        selectionMasks['n_cls_dt']     =ak.num(dt_cluster,axis=1)>=1
        
        selectionMasks["pid_gamma"]    =( clusterMasks.pid_gamma )
        selectionMasks["pid_gluon"]    =( clusterMasks.pid_gluon )
        selectionMasks["pid_quarks"]    =( clusterMasks.pid_quarks )
        #selectionMasks["pid_muon"]    =( clusterMasks.pid_muons )
        selectionMasks["pid_mesons"]    =( clusterMasks.pid_mesons )

        CSC_sel_ABCD = ["ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto",
                        "IntimeCut","timeSpreadCut","ClusterID","genMuon"]
        CSC_sel_OOT  = ["ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto",
                        "OOT_timeCut","timeSpreadCut","ClusterID"]
        CSC_sel_ABCD_inv_ME11 = ["ME11_12_veto_inv","jetVeto","muonVeto","MB1seg_veto","RB1_veto",
                        "IntimeCut","timeSpreadCut","ClusterID","genMuon"]
        if isData:
            CSC_sel_ABCD_inv_ME11 = ["ME11_12_veto_inv","jetVeto","muonVeto","MB1seg_veto","RB1_veto",
                                     "IntimeCut","timeSpreadCut","ClusterID"]
        CSC_sel_ABCD_data = ["ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto",
                        "IntimeCut","timeSpreadCut","ClusterID","NRechit100"]
        

        selectionMasks['cls_ABCD']  = buildMask(clusterMasks,CSC_sel_ABCD)
        selectionMasks['cls_OOT']   = buildMask(clusterMasks,CSC_sel_OOT)
        selectionMasks['cls_ABCD_data']  = buildMask(clusterMasks,CSC_sel_ABCD_data)


        selectionMasks['cls_StatVeto']     =  buildMask(clusterMasks,['ME11_12_veto','MB1seg_veto','RB1_veto'])     
        selectionMasks['cls_JetMuVeto']    =  buildMask(clusterMasks,['jetVeto','muonVeto'])                
        selectionMasks['cls_JetMuStaVeto'] =  buildMask(clusterMasks,['jetVeto','muonVeto','ME11_12_veto','MB1seg_veto','RB1_veto']) 
        selectionMasks['cls_ABCD_inv_ME11']  = buildMask(clusterMasks,CSC_sel_ABCD_inv_ME11)

        DT_sel_OOT  = ["dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_OOT",]
       
        DT_sel_ABCD = ["dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time","dt_genMuon","dt_deadzones"]
      
        DT_sel_ABCD_inv_MB1 = ["dt_MB1veto_inv","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time","dt_deadzones","dt_genMuon"]
     
        if isData:
            DT_sel_ABCD_inv_MB1 = ["dt_MB1veto_inv","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time","dt_deadzones"]
    
        DT_sel_ABCD_data = ["dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time","NRechit70","dt_deadzones"]

       
        DT_sel_ABCD_inv_jet_veto = ["dt_MB1veto","dt_jetVeto_inv","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time"]
        DT_sel_ABCD_old_Muon_Veto = ["dt_jetVeto","MuonVeto_0p4","dt_MB1veto","dt_RPC","dt_MB1adj","dt_time"]
        DT_sel_ABCD_new_Muon_Veto = ["dt_jetVeto","MuonVeto_0p8","dt_MB1veto","dt_RPC","dt_MB1adj","dt_time"]
        #only_old_Muon_Veto = ["MuonVeto_0p4"]
        #only_new_Muon_Veto = ["MuonVeto_0p8"]

        selectionMasks['dt_cls_OOT']  = buildMask(dt_clusterMasks,DT_sel_OOT)         
        selectionMasks['dt_cls_ABCD']  = buildMask(dt_clusterMasks,DT_sel_ABCD)         
        selectionMasks['dt_cls_ABCD_data']  = buildMask(dt_clusterMasks,DT_sel_ABCD_data)
        

       
        #selectionMasks['dt_cls_ABCD_old_Muon_Veto']  = buildMask(dt_clusterMasks,DT_sel_ABCD_old_Muon_Veto)
        #selectionMasks['dt_cls_ABCD_new_Muon_Veto']  = buildMask(dt_clusterMasks,DT_sel_ABCD_new_Muon_Veto)

      
        selectionMasks['dt_cls_ABCD_inv_MB1']  = buildMask(dt_clusterMasks,DT_sel_ABCD_inv_MB1)

        selectionMasks['dt_cls_ABCD_inv_jet_veto']  = buildMask(dt_clusterMasks,DT_sel_ABCD_inv_jet_veto)
        
        #selectionMasks["dt_only_old_Muon_Veto"] = buildMask(dt_clusterMasks, only_old_Muon_Veto )
        #selectionMasks["dt_only_new_Muon_Veto"] = buildMask(dt_clusterMasks, only_new_Muon_Veto )

        if self.isElectronChannel:
            preselections = ['trigger_ele','MET',"METfilters",'good_lepton']       
        else:
            preselections = ['trigger_mu',
                             'MET',
                             "METfilters",
                             'good_lepton'
                             ]       
            preselections_no_met = ['trigger_mu',
                                      'good_lepton'  
                             ]
        regions = {
            "PreSel_dt"       :preselections,            
            "PreSel"          :preselections,
            #"ele_W_CR"     :['trigger_ele','MET',"METfilters",'good_electron',"W_CR",],
           # "JetMuVeto"    :preselections+["cls_JetMuVeto"],
           # "JetMuStaVeto" :preselections+["cls_JetMuStaVeto"],
            #"ABCD"         :preselections+["cls_ABCD"],            
           #"ABCD_old_Muon_Veto"         :preselections+["cls_ABCD_old_Muon_Veto"],
            #"ABCD_OOT"     :preselections+["cls_OOT"],
            #"ABCD_dt"      :preselections+["dt_cls_ABCD"],            
            #"ABCD_dt_old_Muon_Veto"      :preselections+["dt_cls_ABCD_old_Muon_Veto"],
            #"ABCD_dt_new_Muon_Veto"      :preselections+["dt_cls_ABCD_new_Muon_Veto"],
            #"dt_only_old_Muon_Veto"      : ["dt_only_old_Muon_Veto"],
            #"dt_only_new_Muon_Veto"      : ["dt_only_new_Muon_Veto"],
            "ABCD_dt_OOT"  :preselections+["dt_cls_OOT"],
            "ABCD_dt_MC" : preselections+["dt_cls_ABCD"],     
            "ABCD_dt_inv_MB1" : preselections+["dt_cls_ABCD_inv_MB1"], 
            "ABCD_dt_Data": preselections+["dt_cls_ABCD_data"],            
            
           
            "ABCD_OOT"  :preselections+["cls_OOT"],
            "ABCD_MC" : preselections+["cls_ABCD"], 
            "ABCD_inv_ME11" : preselections+["cls_ABCD_inv_ME11"], 
            "ABCD_data" : preselections+["cls_ABCD_data"], 
          
            "ABCD_dt_OOT_no_met"  :preselections_no_met+["dt_cls_OOT"],
            "ABCD_dt_MC_no_met" : preselections_no_met+["dt_cls_ABCD"],
            "ABCD_dt_inv_MB1_no_met" : preselections_no_met+["dt_cls_ABCD_inv_MB1"],
            "ABCD_dt_Data_no_met" : preselections_no_met+["dt_cls_ABCD_data"],
            
            "ABCD_OOT_no_met"  :preselections_no_met+["cls_OOT"],
            "ABCD_MC_no_met" : preselections_no_met+["cls_ABCD"],   
            "ABCD_inv_ME11_no_met" : preselections_no_met+["cls_ABCD_inv_ME11"],
            "ABCD_Data_no_met" : preselections_no_met+["cls_ABCD_data"],  



           #"ABCD_dt_inv_MB1_pick_events" : preselections+["dt_cls_ABCD_inv_MB1_pick_events"], 
           #"ABCD_inv_ME11"   : preselections+["cls_ABCD_inv_ME11"],
            # "ABCD_dt_inv_jet_veto" : preselections+["dt_cls_ABCD_inv_jet_veto"],
        
          # "JetMuStaVeto_dt" :preselections+["dt_JetMuStaVeto"],
            ##"1cls"         :preselections+["n_cls"],            
            #"StatVeto"     :preselections+["cls_StatVeto"],
            #"noselection":[],
        }

        weights = Weights(len(events))
        if not isData: 
            corrections.add_pileup_weight(weights, events.npu,'2018')
            corrections.add_Wpt_kfactor(weights, events.gWPt, dataset)
            if isSignal and "rwctau" in dataset:
                ## expect dataset = "HNL_*_pl{ctau_old}_rwctau{ctau_new}"
                ctau_old =  float(dataset.split("_")[-2].replace("pl",""))/10     ##ctau in dataset name is in mm
                ctau_new =  float(dataset.split("_")[-1].replace("rwctau",""))/10 ##new ctau needs to be in cm
                corrections.add_ctau_weight(weights, llp.ctau, ctau_old, ctau_new)
            pass

        if self._debug:
            print(dataset)
            print("Weight statistics: %r" % weights.weightStatistics) 

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
            cut = selectionMasks["Acceptance_csc"] 

            #Events with clusterID pass
            #llp_selection = maskAndFill(llp.e,ak.any(cluster[cut].llp_match,axis=1),len(llp.e[0])*[0])
            #Events with any cluster matching to llp
            llp_selection = ak.values_astype( ak.any(cluster.llp_match,axis=1),np.int )

            output['llp_cls_eff_z'].fill(dataset=dataset,selection=llp_selection[cut],z=ak.flatten(abs(llp.z[cut])),weight=weights.weight()[cut])
            output['llp_cls_eff_r'].fill(dataset=dataset,selection=llp_selection[cut],r=ak.flatten(llp.r[cut]),weight=weights.weight()[cut])
            output['llp_cls_eff_e'].fill(dataset=dataset,selection=llp_selection[cut],e=ak.flatten(llp.e[cut]),weight=weights.weight()[cut])

            cut = selectionMasks["Acceptance_dt"] 
            llp_selection = ak.values_astype( ak.any(dt_cluster.llp_match,axis=1),np.int )
            output['llp_cls_dt_eff_z'].fill(dataset=dataset,selection=llp_selection[cut],z=ak.flatten(abs(llp.z[cut])),weight=weights.weight()[cut])
            output['llp_cls_dt_eff_r'].fill(dataset=dataset,selection=llp_selection[cut],r=ak.flatten(llp.r[cut]),weight=weights.weight()[cut])
            output['llp_cls_dt_eff_e'].fill(dataset=dataset,selection=llp_selection[cut],e=ak.flatten(llp.e[cut]),weight=weights.weight()[cut])


        output["metXYCorr"].fill(dataset=dataset,region="noselection",metXYCorr=events.metXYCorr,weight=weights.weight()) 
        gparticle_list = ["gamma", "gluon", "quarks","muon","mesons"] 
        
        
        
        #cut = buildMask(selectionMasks,regions["ABCD_dt"]) 
        #if cut.ndim==1:
        #    ev_cut = cut                  ##This is a per-event cut
        #else:
        #    ev_cut = ak.any(cut,axis=1)   ##This is a per-cluster cut, require at least 1 passing cluster

        #ev_cut = ak.fill_none(ev_cut,False)
        #print(os.getcwd())
        #with open('Pick_events.txt', 'w') as f:
        #    f.write('runNum , evtNum , lumiSec \n')
        #    for i in range(len(events.runNum[ev_cut])):
        #        f.write(str(events.runNum[ev_cut][i])+" , "+ str(events.evtNum[ev_cut][i]) + " , "+str(events.lumiSec[ev_cut][i])+"\n")
        #        print(str(events.runNum[ev_cut][i])+" , "+ str(events.evtNum[ev_cut][i]) + " , "+str(events.lumiSec[ev_cut][i])+"\n")
                #print(cluster[cut].dphi_cluster_lep)
                #print(cluster[cut].size)
        ## Fill regions plot
        for region,cutnames in regions.items():
            if "ABCD_dt_MC" in region or "ABCD_MC" in region:
                if isData:
                    continue
            #print(region)
            ## Fill other regions without cutflows
            cut = buildMask(selectionMasks,cutnames)

            if cut.ndim==1:
                ev_cut = cut                  ##This is a per-event cut
            else:
                ev_cut = ak.any(cut,axis=1)   ##This is a per-cluster cut, require at least 1 passing cluster
 
            ev_cut = ak.fill_none(ev_cut,False)
            w_evt = weights.weight()[ev_cut]
            output['weight_histo'].fill(dataset=dataset ,region=region, weight_input  = w_evt  ) 

            if not "dt" in region:
                w_cls      = (weights.weight() * ak.ones_like(cluster.size))[cut] ## use size to pick-up the cluster shape
                output["dphi_cluster_csc"].fill(dataset=dataset,region=region,
                                                ClusterSize=ak.flatten(cluster[cut].size),
                                                dphi_lep =np.abs(ak.flatten(cluster[cut].dphi_cluster_lep)),
                                                dphi_MET=np.abs(ak.flatten(cluster[cut].dphi_cluster_MET)), 
                                                weight=ak.flatten(w_cls))
                output["ClusterID_csc"].fill(dataset=dataset,region=region,
                                                ClusterEta=ak.flatten(cluster[cut].eta),
                                                ClusterAvgStation10 =np.abs(ak.flatten(cluster[cut].AvgStation10)),
                                                ClusterNStation10=np.abs(ak.flatten(cluster[cut].NStation10)),
                                                weight=ak.flatten(w_cls))
                output["dr_lep_and_dphi_lep"].fill(dataset=dataset,region=region,
                                                      dr_lep=np.abs(ak.flatten(cluster[cut].dr_cluster_lep)),
                                                      dphi_lep =np.abs(ak.flatten(cluster[cut].dphi_cluster_lep)),
                                                      weight=ak.flatten(w_cls))
                output["ClusterSize"].fill(dataset=dataset,region=region,
                                           ClusterSize=ak.flatten(cluster[cut].size),
                                           weight=ak.flatten(w_cls))        
                output["ClusterTime"].fill(dataset=dataset,region=region,
                                           ClusterTime=ak.flatten(cluster[cut].time),
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
                output["metXYCorr"].fill(dataset=dataset,region=region,
                                         metXYCorr=events[ev_cut].metXYCorr,
                                        weight=w_evt)                    
                output['jetPt'].fill(dataset=dataset, region = region, 
                                    jetPt = ak.to_numpy(ak.firsts(events[ev_cut].jetPt)),
                                    weight=w_evt)
                output["MT"].fill(dataset=dataset,region=region,MT=events[ev_cut].MT,weight=w_evt)       
            else:
                w_cls      = (weights.weight() * ak.ones_like(dt_cluster.size))[cut] ## use size to pick-up the cluster shape

                output["dphi_cluster_dt"].fill(dataset=dataset,region=region,
                                                ClusterSize=ak.flatten(dt_cluster[cut].size),
                                                dphi_lep =np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_lep)),
                                                dphi_MET=np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_MET)),
                                               
                                                weight=ak.flatten(w_cls))
                output["ClusterID_dt"].fill(dataset=dataset,region=region,
                                                
                                                ClusterEta=ak.flatten(dt_cluster[cut].eta),
                                                ClusterAvgStation10 =np.abs(ak.flatten(dt_cluster[cut].AvgStation10)),
                                                ClusterNStation10=np.abs(ak.flatten(dt_cluster[cut].NStation10)),
                                                weight=ak.flatten(w_cls))
                output["dr_lep_and_dphi_lep_dt"].fill(dataset=dataset,region=region,
                                           dr_lep=np.abs(ak.flatten(dt_cluster[cut].dr_cluster_lep)),
                                           dphi_lep =np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_lep)),
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

                output["ClusterEta_dt_dphi_lep"].fill(dataset=dataset,region=region,
                                           ClusterSize=ak.flatten(dt_cluster[cut].size),
                                           ClusterEta=np.abs(ak.flatten(dt_cluster[cut].eta)),
                                           dphi_lep =np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_lep)),
                                           weight=ak.flatten(w_cls))
                output["ClusterAvgStation10_dt_dphi_lep"].fill(dataset=dataset,region=region,
                                           ClusterSize=ak.flatten(dt_cluster[cut].size),
                                           ClusterAvgStation10=np.abs(ak.flatten(dt_cluster[cut].AvgStation10)),
                                           dphi_lep =np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_lep)),
                                           weight=ak.flatten(w_cls))
                output["ClusterNStation10_dt_dphi_lep"].fill(dataset=dataset,region=region,
                                           ClusterSize=ak.flatten(dt_cluster[cut].size),
                                           ClusterNStation10=ak.flatten(dt_cluster[cut].NStation10),
                                           dphi_lep =np.abs(ak.flatten(dt_cluster[cut].dphi_cluster_lep)),
                                           weight=ak.flatten(w_cls))
                output["ClusternMB1"].fill(dataset=dataset,region=region,
                                           ClusternMB1=ak.flatten(dt_cluster[cut].nMB1),
                                           weight=ak.flatten(w_cls))
                output["Cluster_MuonVeto_Pt_dt"].fill(dataset=dataset,region=region,
                                           Cluster_MuonVeto_pt=ak.flatten(dt_cluster[cut].MuonVetoPt),
                                           weight=ak.flatten(w_cls))
                output["Cluster_MuonVeto_Eta_dt"].fill(dataset=dataset,region=region,
                                           Cluster_MuonVeto_eta=ak.flatten(dt_cluster[cut].MuonVetoEta),
                                           weight=ak.flatten(w_cls))
                output["Cluster_MuonVeto_Phi_dt"].fill(dataset=dataset,region=region,
                                           Cluster_MuonVeto_phi=ak.flatten(dt_cluster[cut].MuonVetoPhi),
                                           weight=ak.flatten(w_cls))
                output["ClusterMaxStationRatio_dt"].fill(dataset=dataset,region=region,
                                           ClusterMaxStationRatio=ak.flatten(dt_cluster[cut].MaxStationRatio),
                                           weight=ak.flatten(w_cls))
                #output["ClusterJetVetoPt_dt"].fill(dataset=dataset,region=region,
                #                           JetVetoPt=ak.flatten(dt_cluster[cut].JetVetoPt),
                #                           weight=ak.flatten(w_cls))
                #output["ClusterJetVetoE_dt"].fill(dataset=dataset,region=region,
                #                           JetVetoE=ak.flatten(dt_cluster[cut].JetVetoE),
                #                           weight=ak.flatten(w_cls))
                #output["ClusterMaxStation_dt"].fill(dataset=dataset,region=region,
                #                           MaxStation=ak.flatten(dt_cluster[cut].MaxStation),
                #                           weight=ak.flatten(w_cls))
                """output["ClusterMaxStationRatio_dt"].fill(dataset=dataset,region=region,
                                           MaxStationRatio=ak.flatten(dt_cluster[cut].MaxStationRatio),
                                           weight=ak.flatten(w_cls))
                output["ClusterNChamber_dt"].fill(dataset=dataset,region=region,
                                           NChamber=ak.flatten(dt_cluster[cut].NChamber),
                                           weight=ak.flatten(w_cls))
                output["ClusterMaxChamber_dt"].fill(dataset=dataset,region=region,
                                           MaxChamber=ak.flatten(dt_cluster[cut].MaxChamber),
                                           weight=ak.flatten(w_cls))
                output["ClusterrpcBx_dt"].fill(dataset=dataset,region=region,
                                           rpcBx=ak.flatten(dt_cluster[cut].rpcBx),
                                           weight=ak.flatten(w_cls))
                output["ClusternMB1_dt"].fill(dataset=dataset,region=region,
                                           nMB1=ak.flatten(dt_cluster[cut].nMB1),
                                           weight=ak.flatten(w_cls))
                output["ClusternRPC_dt"].fill(dataset=dataset,region=region,
                                           nRPC=ak.flatten(dt_cluster[cut].nRPC),
                                           weight=ak.flatten(w_cls))
                output["ClusternMB1_cosmic_minus_dt"].fill(dataset=dataset,region=region,
                                           nMB1_cosmic_minus=ak.flatten(dt_cluster[cut].nMB1_cosmic_minus),
                                           weight=ak.flatten(w_cls))
                output["ClusternMB1_cosmic_plus_dt"].fill(dataset=dataset,region=region,
                                           nMB1_cosmic_plus=ak.flatten(dt_cluster[cut].nMB1_cosmic_plus),
                                           weight=ak.flatten(w_cls))
                """
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
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow_old_Muon_Veto",cutflow=sel,weight=weights.weight()[allcuts])
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow_new_Muon_Veto",cutflow=sel,weight=weights.weight()[allcuts])
        allPreSel= buildMask(selectionMasks,acc_csc_preselections) 
        allPreSel_dt= buildMask(selectionMasks,acc_dt_preselections) 
        
        #for i,sel in enumerate(CSC_sel_ABCD):
        #    allcuts= ak.any( (buildMask(clusterMasks,CSC_sel_ABCD[0:i+1]) & allPreSel), axis=1)     ## select all cuts up to this cut
        #    output['cutflow'].fill(dataset=dataset,region="csc_cutflow",cutflow=sel,weight=weights.weight()[allcuts])
        #for i,sel in enumerate(DT_sel_ABCD):
        #    allcuts= ak.any( (buildMask(dt_clusterMasks,DT_sel_ABCD[0:i+1]) & allPreSel_dt), axis=1)     ## select all cuts up to this cut
        #    output['cutflow'].fill(dataset=dataset,region="dt_cutflow",cutflow=sel,weight=weights.weight()[allcuts])

        """for i,sel in enumerate(DT_sel_ABCD_old_Muon_Veto):
            allcuts= ak.any( (buildMask(dt_clusterMasks,DT_sel_ABCD_old_Muon_Veto[0:i+1]) & allPreSel_dt), axis=1)     ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow_old_Muon_Veto",cutflow=sel,weight=weights.weight()[allcuts])
        for i,sel in enumerate(DT_sel_ABCD_new_Muon_Veto):
            allcuts= ak.any( (buildMask(dt_clusterMasks,DT_sel_ABCD_new_Muon_Veto[0:i+1]) & allPreSel_dt), axis=1)     ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow_new_Muon_Veto",cutflow=sel,weight=weights.weight()[allcuts])
        for i,sel in enumerate(only_old_Muon_Veto):
            allcuts= ak.any( (buildMask(dt_clusterMasks,only_old_Muon_Veto[0:i+1]) & allPreSel_dt), axis=1)     ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow_only_old_Muon_Veto",cutflow=sel,weight=weights.weight()[allcuts])
        for i,sel in enumerate(only_old_Muon_Veto):
            allcuts= ak.any( (buildMask(dt_clusterMasks,only_old_Muon_Veto[0:i+1]) & allPreSel_dt), axis=1)     ## select all cuts up to this cut
            output['cutflow'].fill(dataset=dataset,region="dt_cutflow_only_old_Muon_Veto",cutflow=sel,weight=weights.weight()[allcuts])
        """
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
            destination = "root://cmseos.fnal.gov//store/user/lpclonglived/HNL/"

            cut = ak.any(
                    buildMask(selectionMasks,regions["ABCD"]) 
                    ,axis=1)
            cut = cut | ak.any(buildMask(selectionMasks,regions["ABCD_dt"]),axis=1)
            if ak.any(cut,axis=0):
                print("Found events pass skim cut, writing out")
                with uproot.recreate(filename) as fout:
                    cluster['passABCD'] =  buildMask(selectionMasks,regions["ABCD"])
                    dt_cluster['passABCD_dt'] =  buildMask(selectionMasks,regions["ABCD_dt"])
                    #fout["MuonSystem"] = uproot_writeable(events[cut], events.fields)    # TODO: find out why we can't write all event fields
                    fout["MuonSystem"] = {"cluster":cluster[cut],"dt_cluster":dt_cluster[cut],"gParticle":gParticle[cut]}
            
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

    def postprocess(self, accumulator):
        # set everything to 1/fb scale
        lumi = 1000  # [1/pb]

        scale = {}
        for dataset, dataset_sumw in accumulator['sumw'].items():
            if "rwctau" in dataset:
                ctau    = float(dataset.split("_")[-1].replace("rwctau",""))
                mass    = dataset.split("_")[2].replace("mHNL","")
                xsec    = corrections.reweightXsec(ctau,mass)
                #print("reweighting ct for ",dataset," with xsec = ",xsec)
                scale[dataset] = lumi*xsec/dataset_sumw
                
            elif dataset in corrections.load_xsection().keys():
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

