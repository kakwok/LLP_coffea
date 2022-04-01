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
        
        pid_gamma     = (cluster.Cluster_match_gParticle_id==22)
        pid_gluon     = (cluster.Cluster_match_gParticle_id==21)
        pid_quarks    = (abs(cluster.Cluster_match_gParticle_id)<10)
        pid_muons     = (abs(cluster.Cluster_match_gParticle_id)==13)
        pid_mesons    = (abs(cluster.Cluster_match_gParticle_id)>100) & (abs(cluster.Cluster_match_gParticle_id)<1000)
       
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
            "dr_lep"        : dr_lep        ,
            "pid_gamma"     : pid_gamma     ,
            "pid_gluon"     : pid_gluon     ,
            "pid_quarks"    : pid_quarks    ,
            "pid_muons"     : pid_muons     ,
            "pid_mesons"    : pid_mesons    ,
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
        dr_dt_cluster_lep = ak.fill_none(dt_cluster_dir.delta_r(ak.firsts(good_lep)),-999)
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
                 "dr_cluster_lep":dr_dt_cluster_lep,
                 "Cluster_match_gParticle_id":events.dtRechitCluster_match_gParticle_Id,
                 "Cluster_match_gParticle_minDeltaR":events.dtRechitCluster_match_gParticle_deltaR, 
                 "Cluster_match_gParticle_eta":events.dtRechitCluster_match_gParticle_Eta,
                 "Cluster_match_gParticle_phi":events.dtRechitCluster_match_gParticle_Phi,
                 "Cluster_match_gParticle_e":events.dtRechitCluster_match_gParticle_E,
                 "Cluster_match_gParticle_pt":events.dtRechitCluster_match_gParticle_Pt,
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
        dt_dphi_MET  = (abs(dt_cluster.dphi_cluster_MET)<0.75)
        dt_size      = (dt_cluster.size>=100)
        dr_lep      = (dt_cluster.dr_cluster_lep>0.4)
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
                "dr_lep"     :dr_lep     ,
        })
        return clusterMasks
    
    def process(self, events):
        output = self.accumulator.identity()  ## get from histograms
        dataset = events.metadata['dataset']        
        
        start,stop = events._branchargs['entry_start'],events._branchargs['entry_stop']
        events = uproot.lazy(events._tree)
        events = events[start:stop]
                
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
                                
                  
       ## All possible pairs 
        #cls_lep_pair = ak.cartesian({"cls":cluster_dir,'lep':lep},axis=1,nested=True)
        #dphi_lep_cls = cls_lep_pair.cls.delta_phi(cls_lep_pair.lep)       

        llp      = self.buildLLP(events)
        good_lep,ele,muons = self.buildGoodLeptons(events) 
        cluster = self.buildCSCcluster(events,good_lep)        
        dt_cluster = self.buildDTcluster(events,good_lep)        

        clusterMasks = self.selectCSCcluster(cluster,events) 
        dt_clusterMasks = self.selectDTcluster(dt_cluster,events) 

        #dictionary of cutName:masks
        selectionMasks =   {}

        selectionMasks['Acceptance']   =ak.firsts(events.gLLP_csc)==1
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
        selectionMasks["pid_muon"]    =( clusterMasks.pid_muons )
        selectionMasks["pid_mesons"]    =( clusterMasks.pid_mesons )

        selectionMasks['cls_ABCD']  = (   ((clusterMasks.jetVeto_mask) &(clusterMasks.muonVeto_mask))
            & (clusterMasks.ME11_12_veto)
            &((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto))
            & (clusterMasks.IntimeCut)
            & (clusterMasks.timeSpreadCut)
            & (clusterMasks.ClusterID))
        selectionMasks["cls_OOT"] =       (     ( clusterMasks.jetVeto_mask &clusterMasks.muonVeto_mask)
            & (clusterMasks.ME11_12_veto)
            &((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto))
            & (clusterMasks.OOT_timeCut)
            & (clusterMasks.timeSpreadCut)
            & (clusterMasks.ClusterID))

        selectionMasks['cls_StatVeto']     =  (clusterMasks.ME11_12_veto)& ((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto))     
        selectionMasks['cls_JetMuVeto']    =  (clusterMasks.jetVeto_mask) &(clusterMasks.muonVeto_mask)                
        selectionMasks['cls_JetMuStaVeto'] =( (clusterMasks.jetVeto_mask) &(clusterMasks.muonVeto_mask)
                                            & (clusterMasks.ME11_12_veto)
                                            &((clusterMasks.MB1seg_veto) & (clusterMasks.RB1_veto)))

        selectionMasks['dt_cls_OOT']= (     (( dt_clusterMasks.dt_jetVeto) &(dt_clusterMasks.dt_muonVeto))
                                            & (dt_clusterMasks.dt_MB1veto)
                                            & (dt_clusterMasks.dt_RPC)
                                            & (dt_clusterMasks.dt_MB1adj)
                                            & (dt_clusterMasks.dt_OOT))
        
        selectionMasks['dt_cls_ABCD'] =( (( dt_clusterMasks.dt_jetVeto) &(dt_clusterMasks.dt_muonVeto))
                                            & (dt_clusterMasks.dt_MB1veto)
                                            & (dt_clusterMasks.dt_RPC)
                                            & (dt_clusterMasks.dt_MB1adj)
                                            & (dt_clusterMasks.dt_time))
        if self.isElectronChannel:
            preselections = ['trigger_ele','MET',"METfilters",'good_lepton']       
        else:
            preselections = ['trigger_mu','MET',"METfilters",'good_lepton']       


     
        regions = {
            "PreSel"       :preselections,            
            #"ele_W_CR"     :['trigger_ele','MET',"METfilters",'good_electron',"W_CR",],
            "ABCD"         :preselections+["cls_ABCD"],            
            "ABCD_OOT"     :preselections+["cls_OOT"],
            "ABCD_dt"      :preselections+["dt_cls_ABCD"],            
            "ABCD_dt_OOT"  :preselections+["dt_cls_OOT"],
            ##"1cls"         :preselections+["n_cls"],            
            #"JetMuVeto"    :preselections+["cls_JetMuVeto"],
            #"JetMuStaVeto" :preselections+["cls_JetMuStaVeto"],
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
        
        
        def buildMask(allMasks,cutnames):
            allcuts = allMasks[cutnames[0]]
            for i,cutname in enumerate(cutnames):
                allcuts = (allcuts) & allMasks[cutname]
            return allcuts
        
        
        if isSignal:
            output['accept'].fill(dataset=dataset,
                                  gLLP_csc=ak.firsts(events.gLLP_csc),
                                  gLLP_dt=gLLP_dt,weight=weights.weight()) ## only 1 LLP

            cut = selectionMasks["Acceptance"]
            output['gLLP_e'].fill(dataset=dataset,gLLP_e = ak.firsts(llp[cut].e) , weight=weights.weight()[cut])
            output['gLLP_pt'].fill(dataset=dataset,gLLP_pt = ak.firsts(llp[cut].pt), weight=weights.weight()[cut])
            output['gLLP_eta'].fill(dataset=dataset,gLLP_eta = ak.firsts(llp[cut].eta), weight=weights.weight()[cut])
            output['glepdPhi'].fill(dataset=dataset,gLLP_lepdPhi = np.abs(ak.flatten(events[cut].gLLP_lepdPhi)), weight=weights.weight()[cut])
           
        gparticle_list = ["gamma", "gluon", "quarks","muon","mesons"]
        ## Fill regions plot
        for region,cutnames in regions.items():
        
            ## Fill other regions without cutflows
            cut = buildMask(selectionMasks,cutnames)

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
                output["ClusterID_csc"].fill(dataset=dataset,region=region,
                                                ClusterEta=ak.flatten(cluster[cut].eta),
                                                ClusterAvgStation10 =np.abs(ak.flatten(cluster[cut].AvgStation10)),
                                                ClusterNStation10=np.abs(ak.flatten(cluster[cut].NStation10)),
                                                weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_id"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_id=ak.flatten(cluster[cut].Cluster_match_gParticle_id),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_minDeltaR"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_minDeltaR=ak.flatten(cluster[cut].Cluster_match_gParticle_minDeltaR),
                                           weight=ak.flatten(w_cls)) 
                output["Cluster_match_gParticle_eta"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_eta=ak.flatten(cluster[cut].Cluster_match_gParticle_eta),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_phi"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_phi=ak.flatten(cluster[cut].Cluster_match_gParticle_phi),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_pt"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_pt=ak.flatten(cluster[cut].Cluster_match_gParticle_pt),
                                           weight=ak.flatten(w_cls))                
                output["Cluster_match_gParticle_e"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_e=ak.flatten(cluster[cut].Cluster_match_gParticle_e),
                                           weight=ak.flatten(w_cls))
                for particle_name in gparticle_list:
                    #Add pid req
                    cut_pid = buildMask(selectionMasks,cutnames+["pid_"+str(particle_name)])

                    if cut_pid.ndim==1:
                        ev_cut_pid = cut_pid                  ##This is a per-event cut
                    else:
                        ev_cut_pid = ak.any(cut_pid,axis=1)   ##This is a per-cluster cut, require at least 1 passing cluster 
                   
                    w_cls_pid      = (weights.weight() * ak.ones_like(cluster.size))[cut_pid]
                    
                    output["Cluster_match_gParticle_minDeltaR_cat"].fill(dataset=dataset,region=region,
                                               pdgId=particle_name,
                                               Cluster_match_gParticle_minDeltaR=ak.flatten(cluster[cut_pid].Cluster_match_gParticle_minDeltaR),
                                               weight=ak.flatten(w_cls_pid))
                    output["Cluster_match_gParticle_eta_cat"].fill(dataset=dataset,region=region,
                                               pdgId=particle_name,
                                               Cluster_match_gParticle_eta=ak.flatten(cluster[cut_pid].Cluster_match_gParticle_eta),
                                               weight=ak.flatten(w_cls_pid))
                    output["Cluster_match_gParticle_phi_cat"].fill(dataset=dataset,region=region,
                                               pdgId=particle_name,
                                               Cluster_match_gParticle_phi=ak.flatten(cluster[cut_pid].Cluster_match_gParticle_phi),
                                               weight=ak.flatten(w_cls_pid))
                    output["Cluster_match_gParticle_pt_cat"].fill(dataset=dataset,region=region,
                                               pdgId=particle_name,
                                               Cluster_match_gParticle_pt=ak.flatten(cluster[cut_pid].Cluster_match_gParticle_pt),
                                               weight=ak.flatten(w_cls_pid))
                    output["Cluster_match_gParticle_e_cat"].fill(dataset=dataset,region=region,
                                               pdgId=particle_name,
                                               Cluster_match_gParticle_e=ak.flatten(cluster[cut_pid].Cluster_match_gParticle_e),
                                               weight=ak.flatten(w_cls_pid))

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
                output["Cluster_match_gParticle_id_dt"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_id=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_id),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_minDeltaR_dt"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_minDeltaR=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_minDeltaR),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_eta_dt"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_eta=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_eta),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_phi_dt"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_phi=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_phi),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_pt_dt"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_pt=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_pt),
                                           weight=ak.flatten(w_cls))
                output["Cluster_match_gParticle_e_dt"].fill(dataset=dataset,region=region,
                                           Cluster_match_gParticle_e=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_e),
                                           weight=ak.flatten(w_cls))
                """
                output["Cluster_gParticle_id_dt"].fill(dataset=dataset,region=region,
                                           pdgId=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_id),
                                           Cluster_match_gParticle_minDeltaR=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_minDeltaR),
                                           Cluster_match_gParticle_eta=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_eta),
                                           Cluster_match_gParticle_phi=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_phi),
                                           Cluster_match_gParticle_pt=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_pt),
                                           Cluster_match_gParticle_e=ak.flatten(dt_cluster[cut].Cluster_match_gParticle_e),
                                           weight=ak.flatten(w_cls))
                """
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

