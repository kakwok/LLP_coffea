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
from  HNLprocessor.HNLproc_4 import MyProcessor,buildMask,delta_r_pairs
from  HNLprocessor.histograms import *
import warnings
import pickle
import glob
import os
import XRootD
import XRootD.client


ak.behavior.update(candidate.behavior)


class ttbarProcessor(MyProcessor):
    def __init__(self,isElectronChannel=True,saveSkim=False,debug=False):
        super().__init__(isElectronChannel,saveSkim,debug)
        hdict={
             "nbjet": hist.Hist("Events",dataset_axis,region_axis,
                hist.Bin("nbjet", "nbjet", 4, 0, 4)
            ),     
        }
        histograms.update(hdict)
        self._accumulator = processor.dict_accumulator( histograms )


    def buildbjets(self,events,ele,muon):
        jets = ak.zip(
                {k.replace("jet",""):getattr(events,k) for k in events.fields if k.startswith("jet")}
                ,with_name="PtEtaPhiMLorentzVector", 
                behavior=vector.behavior
               )
        bjets = jets[jets.CISV>0.8484] ## medium bjets
        bjets = ak.with_field(bjets,bjets.Phi,"phi")
        bjets = ak.with_field(bjets,bjets.Eta,"eta")

        print("bjets before",len(bjets))
        dphi_e = ak.fill_none(bjets.delta_phi(ak.firsts(ele)),-999,axis=None)
        dr_e   = ak.fill_none(bjets.delta_phi(ak.firsts(ele)),-999,axis=None)
        dphi_mu = ak.fill_none(bjets.delta_phi(ak.firsts(muon)),-999,axis=None)
        dr_mu   = ak.fill_none(bjets.delta_phi(ak.firsts(muon)),-999,axis=None)

        bjets= ak.with_field(bjets,dphi_e,"dphi_e")
        bjets= ak.with_field(bjets,dr_e,"dr_e")
        bjets= ak.with_field(bjets,dphi_mu,"dphi_mu")
        bjets= ak.with_field(bjets,dr_mu,"dr_mu")

        print("bjets after",len(bjets))
        print("bjets",ak.num(bjets,axis=1))
        return bjets
        

    def process(self, events):
        output = self.accumulator.identity()  ## get from histograms
        dataset = events.metadata['dataset']        
                
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

        ## All possible pairs 
        #cls_lep_pair = ak.cartesian({"cls":cluster_dir,'lep':lep},axis=1,nested=True)
        #dphi_lep_cls = cls_lep_pair.cls.delta_phi(cls_lep_pair.lep)       

        llp                = self.buildLLP(events)
        good_lep,ele,muons = self.buildGoodLeptons(events) 
        bjets              = self.buildbjets(events,good_lep,muons) 

        cluster = self.buildCSCcluster(events,good_lep)        
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
        selectionMasks['n_bjets']      =ak.num(bjets,axis=1)==2

        clusterMasks["neg_ME11_12_veto"] = ~clusterMasks['ME11_12_veto']  #make veto mask
        CSC_sel_ABCD = ["ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto", "IntimeCut","timeSpreadCut","ClusterID"]
        CSC_sel_OOT  = ["ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto", "OOT_timeCut","timeSpreadCut","ClusterID"]
        CSC_sel_negME11 = ["neg_ME11_12_veto","jetVeto","muonVeto","MB1seg_veto","RB1_veto", "IntimeCut","timeSpreadCut","ClusterID"]

        selectionMasks['cls_ABCD']  = buildMask(clusterMasks,CSC_sel_ABCD)
        selectionMasks['cls_OOT']   = buildMask(clusterMasks,CSC_sel_OOT)
        selectionMasks['cls_negME11']   = buildMask(clusterMasks,CSC_sel_negME11)

        dt_clusterMasks["neg_dt_MB1veto"] = ~dt_clusterMasks["dt_MB1veto"] #make veto dt mask
        DT_sel_OOT  = ["dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_OOT"]
        DT_sel_ABCD = ["dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time"]
        DT_sel_negMB1 = ["neg_dt_MB1veto","dt_jetVeto","dt_muonVeto","dt_RPC","dt_MB1adj","dt_time"]

        selectionMasks['dt_cls_OOT']  = buildMask(dt_clusterMasks,DT_sel_OOT)         
        selectionMasks['dt_cls_ABCD'] = buildMask(dt_clusterMasks,DT_sel_ABCD)         
        selectionMasks['dt_cls_negMB1'] = buildMask(dt_clusterMasks,DT_sel_negMB1)         

        #bjets selection
        selectionMasks['bjet_cls_dr']    =  ak.all(
                                                ak.all(delta_r_pairs({"cls":cluster,"bjets":bjets})[0]>0.4,axis=2) ## all pairs dr(bjets,cls)>0.4
                                            ,axis=1) # apply to all bjets in the events 
        selectionMasks['bjet_dt_cls_dr'] =  ak.all(
                                                ak.all(delta_r_pairs({"cls":dt_cluster,"bjets":bjets})[0]>0.4,axis=2) ## all dr(bjets,cls)>0.4
                                            ,axis=1)# apply to all bjets in the events 

        if self.isElectronChannel:
            #preselections = ['trigger_ele','MET',"METfilters",'good_lepton',"n_bjets","bjet_cls_dr","bjet_dt_cls_dr"]
            preselections = ['trigger_ele','MET',"METfilters",'good_lepton',"n_bjets","bjet_cls_dr","bjet_dt_cls_dr"]
        else:
            preselections = ['trigger_mu','MET',"METfilters",'good_lepton',"n_bjets","bjet_cls_dr","bjet_dt_cls_dr"]

        regions = {
            "PreSel"       :preselections,            
            "ABCD"         :preselections+["cls_ABCD"],            
            "ABCD_dt"      :preselections+["dt_cls_ABCD"],            
            "ABCD_negME11"      :preselections+["cls_negME11"],            
            "ABCD_negMB1"       :preselections+["dt_cls_negMB1"],            
        }

        #weights = Weights(len(events))
        weights = Weights(len(events))
        print("length of w",len(weights.weight()))
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

 
        ## Fill regions plot
        for region,cutnames in regions.items():

            ## Fill other regions without cutflows
            if self._debug: 
                print(region,cutnames)
                print(selectionMasks)
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
            destination = "root://cmseos.fnal.gov//store/user/kkwok/llp/HNL/skim/"

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

