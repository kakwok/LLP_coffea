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
from  HNLprocessor.histograms import histograms 
import warnings
import pickle
import glob

ak.behavior.update(candidate.behavior)


def maskAndFill(denom,selection,value):
    numer = ak.mask(denom,selection)
    numer = ak.fill_none(numer, value) #fill none with same structure
    return ak.flatten(numer)

class MyProcessor(processor.ProcessorABC):
    def __init__(self,isElectronChannel=True,debug=False, ):
        self._debug = debug
        self.isElectronChannel = isElectronChannel
        self.isMuonChannel = not(isElectronChannel)
        ##define histograms 
        histograms['sumw']= processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator( histograms )

    @property
    def accumulator(self):
        return self._accumulator

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
        dphi_cluster_lep = ak.fill_none(cluster_dir.delta_phi(ak.firsts(good_lep)),-999)
        dr_cluster_lep = ak.fill_none(cluster_dir.delta_r(ak.firsts(good_lep)),-999)
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
                "dphi_cluster_MET":events.cscRechitCluster3MetEENoiseXYCorr_dPhi,                
                "dphi_cluster_lep":dphi_cluster_lep,                
                "dr_cluster_lep":dr_cluster_lep,                
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
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==4) & (abs(cluster.eta)<1.8))|\
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==3) & (abs(cluster.eta)<1.6))|\
        ((cluster.NStation10==1) &(abs(cluster.AvgStation10)==2) & (abs(cluster.eta)<1.6))
        
        muonVeto_mask = ~((muVeto.pt>20) & abs(muVeto.eta<2.4))
        jetVeto_mask  = ~((jetVeto.pt>10)& abs(jetVeto.eta<2.4))
        RE12_veto     = (cluster.RE12==0)
        MB1seg_veto   = (cluster.MB1seg==0)
        RB1_veto      = (cluster.RB1==0)
        ME11_12_veto  = (cluster.ME11_12==0)

        oot_timecut   = (cluster.time < -12.5) #OOT for data
        IntimeCut     = (cluster.time < 12.5) & (cluster.time>-5) ## In-time otherwise        
        timeSpreadCut = (cluster.timeSpread<20)        
        dphi_met      = (abs(cluster.dphi_cluster_MET)<0.75)        
        dphi_lep      = (abs(cluster.dphi_cluster_lep)>2.5)      
        dr_lep      = (cluster.dr_cluster_lep>0.4)

        clusterMasks = ak.zip({
            "ClusterID"     : ClusterID     ,  
            "muonVeto_mask" : muonVeto_mask ,
            "jetVeto_mask"  : jetVeto_mask  ,
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
            })
        return clusterMasks 


    def buildDTcluster(self,events,good_lep):
        dt_cluster_dir= ak.zip(
        {
                'pt':ak.ones_like(events.dtRechitClusterEta),
                "eta":events.dtRechitClusterEta,
                "phi":events.dtRechitClusterPhi,
                'mass':ak.zeros_like(events.dtRechitClusterEta)
            },with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )
        dphi_dt_cluster_lep = ak.fill_none(dt_cluster_dir.delta_phi(ak.firsts(good_lep)),-999)
        #dphi_dt_cluster_mu = ak.fill_none(dt_cluster_dir.delta_phi(ak.firsts(muons)),-999)
        dt_cluster = ak.zip(
            {
                 "size":events.dtRechitClusterSize,
                 "x":events.dtRechitClusterX,
                 "y":events.dtRechitClusterY,
                 "z":events.dtRechitClusterZ,
                 "eta":events.dtRechitClusterEta,
                 "phi":events.dtRechitClusterPhi,
                 "JetVetoPt":events.dtRechitClusterJetVetoPt,
                 "MuonVetoPt":events.dtRechitClusterMuonVetoPt,
                 "NStation10":events.dtRechitClusterNStation10,
                 "AvgStation10":events.dtRechitClusterAvgStation10,
                 "MaxStation":events.dtRechitClusterMaxStation,
                 "nRPC":events.dtRechitCluster_match_RPChits_dPhi0p5,
                 "nMB1_cosmic_minus":events.dtRechitCluster_match_MB1hits_cosmics_minus,
                 "nMB1_cosmic_plus":events.dtRechitCluster_match_MB1hits_cosmics_plus,
                 "nMB1":events.dtRechitCluster_match_MB1hits_0p5,
                 "rpcBx":events.dtRechitCluster_match_RPCBx_dPhi0p5,
                 "dphi_cluster_MET":events.dtRechitClusterMetEENoise_dPhi,
                 "dphi_cluster_lep":dphi_dt_cluster_lep,
                 #"dphi_cluster_mu":dphi_cluster_mu,
            }
        )
        return dt_cluster
 
    def selectDTcluster(self,dt_cluster,events):
        dt_jetVeto  = (dt_cluster.JetVetoPt<20.0)
        dt_muonVeto = (dt_cluster.MuonVetoPt<10.0)
        dt_MB1veto  = (dt_cluster.nMB1<=1)
        dt_RPC      = (dt_cluster.nRPC>=1)
        dt_MB1adj   = (dt_cluster.nMB1_cosmic_minus<=8) & (dt_cluster.nMB1_cosmic_plus<=8)
        dt_time     = (dt_cluster.rpcBx==0)
        dt_OOT      = (dt_cluster.rpcBx>=-100)&(dt_cluster.rpcBx<0)
        dt_dphi_MET  = (abs(dt_cluster.dphi_cluster_MET)<1)
        dt_size      = (dt_cluster.size>=100)
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
        })
        return clusterMasks
    
    def process(self, events):
        output = self.accumulator.identity()
        dataset = events.metadata['dataset']        
        
        start,stop = events._branchargs['entry_start'],events._branchargs['entry_stop']
        events = uproot.lazy(events._tree)
        events = events[start:stop]
                
        isSignal= ak.any(events.gLLP_csc)        
        isData = not(ak.any(events.gLLP_e) or ak.any(events.gLepE))
        isElectronChannel = False
        isMuonChannel = not(isElectronChannel)
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
                                
        llp=ak.zip({
            'pt':events.gLLP_pt,
            #'EMfrac':events.gLLP_EMFracE,
            'e':events.gLLP_e,
            'eta':events.gLLP_eta,
            'z':events.gLLP_decay_vertex_z ,
            'r':events.gLLP_decay_vertex_r,
            'ctau':events.gLLP_ctau,
        })  
                   
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

        ## All possible pairs 
        #cls_lep_pair = ak.cartesian({"cls":cluster_dir,'lep':lep},axis=1,nested=True)
        #dphi_lep_cls = cls_lep_pair.cls.delta_phi(cls_lep_pair.lep)       
 
        cluster = self.buildCSCcluster(events,good_lep)        
        dt_cluster = self.buildDTcluster(events,good_lep)        

        clusterMasks = self.selectCSCcluster(cluster,events) 
        dt_clusterMasks = self.selectDTcluster(dt_cluster,events) 
        ## DT cluster selections
        selection = PackedSelection(np.uint64)        
        
        selection.add('Acceptance',ak.firsts(events.gLLP_csc)==1)
        selection.add('Acceptance_dt',ak.firsts(events.gLLP_dt)==1)
        selection.add('METfilters',events.Flag2_all==True)
        selection.add('trigger_ele',events.SingleEleTrigger==True)
        selection.add('trigger_mu',events.SingleMuonTrigger==True)
        selection.add('good_electron',ak.num(good_ele,axis=1)==1)
        selection.add('good_mu',ak.num(good_mu,axis=1)==1)
        selection.add('MET',events.metEENoise>=30)
        selection.add('W_CR', (events.MT>70) & (events.MT<90) &(events.metXYCorr>60))
       
        cutflow2 = [
            {"name":"cf2_dr_lep"    ,  "cut":clusterMasks.dr_lep},
            {"name":"cf2_ME11_12"    ,  "cut":clusterMasks.ME11_12_veto},
            {"name":"cf2_JetVeto"    ,  "cut":clusterMasks.jetVeto_mask},
            {"name":"cf2_MuVeto"     ,  "cut":clusterMasks.muonVeto_mask},
            {"name":"cf2_MB1seg"     ,  "cut":clusterMasks.MB1seg_veto},
            {"name":"cf2_RB1"        ,  "cut":clusterMasks.RB1_veto},
            {"name":"cf2_Time"       ,  "cut":clusterMasks.IntimeCut},
            {"name":"cf2_TimeSpread" ,  "cut":clusterMasks.timeSpreadCut},
            {"name":"cf2_ClusID"     ,  "cut":clusterMasks.ClusterID},
            {"name":"cf2_dphi_MET"   ,  "cut":clusterMasks.dphi_MET},
            {"name":"cf2_dphi_lep"   ,  "cut":clusterMasks.dphi_lep},
        ]
        allcuts = (clusterMasks.dr_lep)
        for i,s in enumerate(cutflow2):
            allcuts = (allcuts) & (s["cut"])
            selection.add(s["name"],ak.num(cluster[allcuts],axis=1)>0)
            selection.add(s["name"]+"_nminus1",ak.num(cluster[s['cut']],axis=1)>0)

        dt_cutflow = [
            {"name":"dt_JetVeto"    ,  "cut":dt_clusterMasks.dt_jetVeto},
            {"name":"dt_MuVeto"     ,  "cut":dt_clusterMasks.dt_muonVeto},
            {"name":"dt_MB1veto"    ,  "cut":dt_clusterMasks.dt_MB1veto},
            {"name":"dt_RPC"        ,  "cut":dt_clusterMasks.dt_RPC},
            {"name":"dt_MB1adj"     ,  "cut":dt_clusterMasks.dt_MB1adj},
            {"name":"dt_dphi_MET"   ,  "cut":dt_clusterMasks.dt_dphi_MET},
            {"name":"dt_Size"       ,  "cut":dt_clusterMasks.dt_size},
        ]
        allcuts = (dt_clusterMasks.dt_jetVeto)
        for i,s in enumerate(dt_cutflow):
            allcuts = (allcuts) & (s["cut"])
            selection.add(s["name"],ak.num(dt_cluster[allcuts],axis=1)>0)

        selection.add('n_cls',ak.num(cluster,axis=1)>0)
        selection.add('n_cls_dt',ak.num(dt_cluster,axis=1)>0)
        selection.add('nJet',events.nJets>0)

        dt_cls_OOT = dt_cluster[
            (( dt_clusterMasks.dt_jetVeto) &(dt_clusterMasks.dt_muonVeto))
            & (dt_clusterMasks.dt_MB1veto)
            & (dt_clusterMasks.dt_RPC)
            & (dt_clusterMasks.dt_MB1adj)
            & (dt_clusterMasks.dt_OOT)
        ]
        dt_cls_ABCD = dt_cluster[
            (( dt_clusterMasks.dt_jetVeto) &(dt_clusterMasks.dt_muonVeto))
            & (dt_clusterMasks.dt_MB1veto)
            & (dt_clusterMasks.dt_RPC)
            & (dt_clusterMasks.dt_MB1adj)
            & (dt_clusterMasks.dt_time)
        ]
        cls_OOT = cluster[
            ( (clusterMasks.jetVeto_mask) &(clusterMasks.muonVeto_mask))
            & (clusterMasks.ME11_12_veto)
            &((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto))
            & (clusterMasks.OOT_timeCut)
            & (clusterMasks.timeSpreadCut)
            & (clusterMasks.ClusterID)
        ]
        cls_ABCD = cluster[
             ((clusterMasks.jetVeto_mask) &(clusterMasks.muonVeto_mask))
            & (clusterMasks.ME11_12_veto)
            &((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto))
            & (clusterMasks.IntimeCut)
            & (clusterMasks.timeSpreadCut)
            & (clusterMasks.ClusterID)
        ]

        cls_StaVeto = cluster[
             (clusterMasks.ME11_12_veto)& ((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto))
        ]        
        cls_JetMuVeto = cluster[
            ((clusterMasks.jetVeto_mask) &(clusterMasks.muonVeto_mask))
        ]        
        cls_JetMuStaVeto = cluster[
            ( (clusterMasks.jetVeto_mask) &(clusterMasks.muonVeto_mask))
            & (clusterMasks.ME11_12_veto)
            &((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto))
        ]        

        selection.add('cls_OOT',ak.num(cls_OOT,axis=1)>0)
        selection.add('dt_cls_OOT',ak.num(dt_cls_OOT,axis=1)>0)
        selection.add('cls_StatVeto',ak.num(cls_StaVeto,axis=1)>0)
        selection.add('cls_JetMuVeto',ak.num(cls_JetMuVeto,axis=1)>0)
        selection.add('cls_JetMuStaVeto',ak.num(cls_JetMuStaVeto,axis=1)>0)
        selection.add('cls_ABCD',ak.num(cls_ABCD,axis=1)>0)
        selection.add('dt_cls_ABCD',ak.num(dt_cls_ABCD,axis=1)>0)

        if self.isElectronChannel:
            preselections = ['trigger_ele','MET',"METfilters",'good_electron']       
        elif self.isMuonChannel:
            preselections = ['trigger_mu','MET',"METfilters",'good_mu']      
 
        regions = {
            "PreSel"       :preselections,            
            #"ele_W_CR"     :['trigger_ele','MET',"METfilters",'good_electron',"W_CR",],
            "ABCD"         :preselections+["cls_ABCD"],            
            "ABCD_dt"      :preselections+["dt_cls_ABCD"],            
            "ABCD_OOT"     :preselections+["cls_OOT"],
            "ABCD_dt_OOT"  :preselections+["dt_cls_OOT"],
            "1cls"         :preselections+["n_cls"],            
            "JetMuVeto"    :preselections+["cls_JetMuVeto"],
            "JetMuStaVeto" :preselections+["cls_JetMuStaVeto"],
            "StatVeto"     :preselections+["cls_StatVeto"],
            "ABCD_nminus1" :['trigger_ele','good_electron','MET',"METfilters",'n_cls']+[c['name']+"_nminus1" for c in cutflow2],            
            "noselection":[],
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
        #print(dataset)
        #print("Weight statistics: %r" % weights.weightStatistics) 

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
            cut = selection.all(*["Acceptance"])
            output['gLLP_e'].fill(dataset=dataset,gLLP_e = ak.firsts(llp[cut].e) , weight=weights.weight()[cut])
            output['gLLP_pt'].fill(dataset=dataset,gLLP_pt = ak.firsts(llp[cut].pt), weight=weights.weight()[cut])
            output['gLLP_eta'].fill(dataset=dataset,gLLP_eta = ak.firsts(llp[cut].eta), weight=weights.weight()[cut])
            output['glepdPhi'].fill(dataset=dataset,gLLP_lepdPhi = np.abs(ak.flatten(events[cut].gLLP_lepdPhi)), weight=weights.weight()[cut])

        ## Fill n-1 plot
        cut = selection.all(*set([]))
        nMinus1cuts = [c['name']+"_nminus1" for c in cutflow2]
        if self.isElectronChannel:
            presel = ['trigger_ele','good_electron','MET',"METfilters",'n_cls']
        elif self.isMuonChannel:
            presel = ['trigger_mu','good_mu','MET',"METfilters",'n_cls']

        #output["nCluster_n-1"].fill(dataset=dataset,nCluster=ak.num(cluster[cut],axis=1),Nminus1=0,weight=weights.weight()[cut])
        #cut = selection.all(*presel)
        #output["nCluster_n-1"].fill(dataset=dataset,nCluster=ak.num(cluster[cut],axis=1),Nminus1=1,weight=weights.weight()[cut])
        #for i, cut in enumerate(nMinus1cuts):
        #    nMinus1cut = nMinus1cuts[:i]+nMinus1cuts[i+1:] + presel
        #    cut = selection.all(*nMinus1cut)
        #    output["nCluster_n-1"].fill(dataset=dataset,nCluster=ak.num(cluster[cut],axis=1),Nminus1=i+2,weight=weights.weight()[cut])            
        #cut = selection.all(*regions["ABCD_nminus1"])
        #output["nCluster_n-1"].fill(dataset=dataset,nCluster=ak.num(cluster[cut],axis=1),Nminus1=len(nMinus1cuts)+2,weight=weights.weight()[cut])            

        cf_regions ={
            "signal_ABCD_cf2":["Acceptance"]+preselections+['n_cls']+[c['name'] for c in cutflow2],            
            "ABCD_cf2"       :preselections+['n_cls']+[c['name'] for c in cutflow2],            
            "signal_ABCD_dt" :["Acceptance_dt"]+preselections+['n_cls_dt']+[c['name'] for c in dt_cutflow],            
            "ABCD_dt"        :preselections+['n_cls_dt']+[c['name'] for c in dt_cutflow],            
        }
        elif isElectronChannel:
            cf_regions ={
                "ele_signal_ABCD_cf2":["Acceptance",'trigger_ele','good_electron','MET',"METfilters",'n_cls']+[c['name'] for c in cutflow2_ele],            
                "ele_ABCD_cf2"       :['trigger_ele','good_electron','MET',"METfilters",'n_cls']+[c['name'] for c in cutflow2_ele],            
                "ele_ABCD_cf3"       :['trigger_ele','good_electron','MET',"METfilters",'n_cls']+[c['name'] for c in cutflow3_ele],             
        }

        for region,cuts in cf_regions.items():
            ## Fill cut flow plots 
            allcuts = set([])
            cut = selection.all(*allcuts)
            ### CSC cutflows
            if not "dt" in region:
                try:
                    output["nCluster"].fill(dataset=dataset,region=region,
                                            nCluster=ak.num(cluster[cut],axis=1),
                                            cutFlow=0,
                                            weight=weights.weight()[cut])            
                    ## Fill 1-to-n-th cutflow
                    for i, cut in enumerate(cuts):
                        allcuts.add(cut)
                        cut = selection.all(*allcuts)
                        output["nCluster"].fill(dataset=dataset,region=region,
                                            nCluster=ak.num(cluster[cut],axis=1),
                                            cutFlow=i+1,
                                            weight=weights.weight()[cut])          
                except ValueError:
                    print(ak.num(cluster[cut],axis=1))
            else: 
                ## Fill dt cutflows
                output["nCluster_dt"].fill(dataset=dataset,region=region,
                                        nCluster=ak.num(dt_cluster[cut],axis=1),
                                        cutFlow=0,
                                        weight=weights.weight()[cut])            
                ## Fill 1-to-n-th cutflow
                for i, cut in enumerate(cuts):
                    allcuts.add(cut)
                    cut = selection.all(*allcuts)
                    output["nCluster_dt"].fill(dataset=dataset,region=region,
                                        nCluster=ak.num(dt_cluster[cut],axis=1),
                                        cutFlow=i+1,
                                        weight=weights.weight()[cut])           

        ## Fill regions plot
        for region,cuts in regions.items():

            ## Fill other regions without cutflows
            cut = selection.all(*cuts)
            ##per-event weight
            weight = weights.weight()[cut]

            if not "dt" in region:
                w_cls      = (weights.weight() * ak.ones_like(cluster.size))[cut] ## use size to pick-up the cluster shape
                output["dphi_cluster_lep"].fill(dataset=dataset,region=region,
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
                output["metXYCorr"].fill(dataset=dataset,region=region,
                                         metXYCorr=events[cut].metXYCorr,
                                        weight=weight)                    
                output['jetPt'].fill(dataset=dataset, region = region, 
                                    jetPt = ak.to_numpy(ak.firsts(events[cut].jetPt)),
                                    weight=weight)
                output["MT"].fill(dataset=dataset,region=region,MT=events[cut].MT,weight=weight)       
            else:
                w_cls      = (weights.weight() * ak.ones_like(dt_cluster.size))[cut] ## use size to pick-up the cluster shape

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

 
        return output

    def postprocess(self, accumulator):
        # set everything to 1/fb scale
        lumi = 1000  # [1/pb]

        scale = {}
        for dataset, dataset_sumw in accumulator['sumw'].items():
            if "rwctau" in dataset:
                ctau    = float(dataset.split("_")[-1].replace("rwctau",""))
                mass    = dataset.split("_")[2].replace("mHNL","").replace("p0","")
                xsec    = corrections.reweightXsec(ctau,mass)
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

