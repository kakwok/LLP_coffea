import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import coffea
from coffea import hist, processor
import uproot
from optparse import OptionParser

import pickle
import glob

import mplhep as hep
plt.style.use(hep.style.CMS) 
import pandas as pd

def relabel(h):
    a = h.axis('dataset')
    for s in a.identifiers():
        if "testpoint" in s.name or "=" in s.name or len(s.name.split("_"))==0: continue
        if "HNL" in s.name:
            if "electron" in s.name: flavour = "e"
            elif "muon" in s.name: flavour ="\mu"
            if "rwctau" in s.name:
                m = float(s.name.split("_")[-3].replace("mHNL","").replace("p","."))
                ct = int(s.name.split("_")[-1].replace("rwctau",""))
#                 print("m=%s,ct=%s"%(m,ct)    )                
            else:
                m = float(s.name.split("_")[-2].replace("mHNL","").replace("p","."))
                ct = int(s.name.split("_")[-1].replace("pl",""))
#                 print("m=%s,ct=%s"%(m,ct)    )
            a.index(s).label = r"$m_{HNL}^{%s}=%s GeV,c\tau$=%s m"%(flavour,m,ct/1000)
        elif "WJets" in s.name:
            a.index(s).label = "Wjets"
    return 

def rename(name):
    if "HNL" in name:
        m = name.split("_")[2].replace("mHNL","").replace("p0","GeV")
        if "rwctau" in name:
            ct = (name.split("_")[-1].replace("rwctau",""))
        else:
            ct = (name.split("_")[-1].replace("pl",""))
        label = "_".join([m,ct+"mm"])
    elif "WJet" in name:
        label= "WJetsToLNu"
    elif "EGamma" in name:
        label = "EGamma"
    return label

def loadPickle(f):
    import HNLprocessor.corrections as corrections

    xsections = corrections.load_xsection()

    # with open('../HNL_histograms_Mar15_muons_signal_Wpt.pickle','rb') as f:                
    with open(f,'rb') as f:                
    # with open('../test.pickle','rb') as f:                    
        out = pickle.load(f)
    
    lumi = 137 ## fb
    
    for k,h in out.items():
    #     print(k,)
        if (type(h)!=hist.Hist): continue
        h.scale({ d: lumi for d in h.identifiers("dataset") if not "2018" in d.name}, axis="dataset")
    return out

import matplotlib.patches as patches

def plotEff(h_list,ax,axis='e'):
    for h in h_list:
        h_pass=h.integrate('selection',slice(1,None))
        h_pass.label='Cluster Efficiency'
        
        hist.plotratio(
            ax=ax,
            num  =h_pass.project(axis),
            denom=h.project(axis),
            error_opts={'marker': '.'},
            unc='num',
            label=h.identifiers('dataset')[0].label,
            clear=False
        )

    ax.legend(loc='best')
    return ax
def drawCSCsteel(ax,hORv='v'): 
    y_max = ax.get_ylim()[1]
    if hORv=='v':
        ax.axvline(632);ax.axvline(632+39)
        ax.axvline(724);ax.axvline(724+64)
        ax.axvline(849);ax.axvline(849+62)
        ax.axvline(970);ax.axvline(970+32)
    else:
        ax.axhline(632);ax.axhline(632+39)
        ax.axhline(724);ax.axhline(724+64)
        ax.axhline(849);ax.axhline(849+62)
        ax.axhline(970);ax.axhline(970+32)   
        
    ax.text(570, y_max*1.02, 'ME1/1', fontsize=12)
    ax.text(670, y_max*1.02, 'ME1/2-3', fontsize=12)
    ax.text(800, y_max*1.02, 'ME2', fontsize=12)
    ax.text(920, y_max*1.02, 'ME3', fontsize=12)
    ax.text(1015, y_max*1.02,'ME4', fontsize=12)        
    return 

def drawCSCr(ax):
    y_max = ax.get_ylim()[1]
    ax.axvline(350,linestyle="--",color='grey')
    ax.text(350-110,y_max*0.05, "Inner ring", fontsize=15)
    ax.text(350+15 ,y_max*0.05, "Outer ring", fontsize=15)
    return ax

def drawCSCz(ax,text_loc=0.7):    
    ax.set_xlim(400,1100)
    (xmin,xmax) = ax.get_xlim()

    y_max = ax.get_ylim()[1]

    preME11 = patches.Rectangle((xmin, 0), 568-xmin, 2,color='grey',alpha=0.3)
    ME11_12 = patches.Rectangle((632, 0), 39, 2,color='grey',alpha=0.3)
    ME12_2  = patches.Rectangle((724, 0), 65, 2,color='grey',alpha=0.3)
    ME2_3   = patches.Rectangle((849, 0), 62, 2,color='grey',alpha=0.3)
    ME3_4   = patches.Rectangle((970, 0), 32, 2,color='grey',alpha=0.3)
    beyond  = patches.Rectangle((1050, 0),50, 2,color='grey',alpha=0.3)

    ax.text(570, y_max*1.02, 'ME1/1', fontsize=12)
    ax.text(670, y_max*1.02, 'ME1/2-3', fontsize=12)
    ax.text(800, y_max*1.02, 'ME2', fontsize=12)
    ax.text(920, y_max*1.02, 'ME3', fontsize=12)
    ax.text(1015, y_max*1.02,'ME4', fontsize=12)
    ax.text(xmin+5 ,y_max*0.15, "Steel", fontsize=15,rotation=90)
    ax.text(xmax-20,y_max*0.15, "Beyond CMS", fontsize=15,rotation=90)

    ax.add_patch(preME11)
    ax.add_patch(ME11_12)
    ax.add_patch(ME12_2)
    ax.add_patch(ME2_3)
    ax.add_patch(ME3_4)
    ax.add_patch(beyond)
    return ax


def drawDTz(ax):
    y_max = ax.get_ylim()[1]
    ax.axvline(140,linestyle="--",color='grey')
    ax.axvline(400,linestyle="--",color='grey')    

    ax.text(140-110,y_max*0.05, "Wheel 0", fontsize=15)
    ax.text(220 ,y_max*0.05, "Wheel 1", fontsize=15)
    ax.text(500 ,y_max*0.05, "Wheel 2", fontsize=15)    
    return ax

def drawDTr(ax,text_loc=0.7):    
    ax.set_xlim(200,800)
    (xmin,xmax) = ax.get_xlim()

    y_max = ax.get_ylim()[1]

    preMB1 = patches.Rectangle((xmin, 0), 405-xmin, 2,color='grey',alpha=0.3)
    MB1_2  = patches.Rectangle((465, 0) , 20, 2,color='grey',alpha=0.3)
    MB2_3  = patches.Rectangle((540, 0) , 55, 2,color='grey',alpha=0.3)
    MB3_4  = patches.Rectangle((640, 0) , 60, 2,color='grey',alpha=0.3)
    

    ax.text(420, y_max*1.02, 'MB1', fontsize=12)
    ax.text(500, y_max*1.02, 'MB2', fontsize=12)
    ax.text(600, y_max*1.02, 'MB3', fontsize=12)
    ax.text(710, y_max*1.02, 'MB4', fontsize=12)
    ax.axvline(740,linestyle="--",color='grey')    
    
    ax.text(xmin+5 ,y_max*0.15, "Steel/Solenoid", fontsize=15,rotation=90)
    ax.text(xmax-10,y_max*0.15, "Beyond CMS", fontsize=15,rotation=90)
    

    ax.add_patch(preMB1)
    ax.add_patch(MB1_2)
    ax.add_patch(MB2_3)
    ax.add_patch(MB3_4)
    return ax

def plotEffrz(hr,hz,datasets,outf,drawCSC="CSC",ylim=None,close=True):
    figsize=(18,6)
    fig, axs = plt.subplots(1,2,figsize=figsize)
    
    axs= axs.flatten()
    relabel(hz)
    relabel(hr)
    
    plotEff([hr[d] for d in datasets],axs[0],"r")
    plotEff([hz[d] for d in datasets],axs[1],"z")
    if drawCSC=="CSC":
        drawCSCz(axs[1])
        drawCSCr(axs[0])
    else:
        drawDTz(axs[1])
        drawDTr(axs[0])
    axs[0].legend(loc=1)
    plt.tight_layout()
    fig.savefig(outf)
    if close:
        plt.close(fig)
    return fig,axs

def plotEff_e(he,datasets,outf,close=True):
    
    figsize=(12,8)
    fig, axs = plt.subplots(1,1,figsize=figsize)
    
    relabel(he)
    
    plotEff([he[d] for d in datasets],axs,"e")
    hep.cms.label(ax=axs)
    fig.savefig(outf)
    if close:
        plt.close(fig)
    return fig,axs

   
def plot(h,dataset,outf,density=True,logy=None,ylim=None,xlim=None,close=True):

    plt.style.use(hep.style.CMS) 
    
    fig, axs = plt.subplots(1,1, figsize=(12,8))
    relabel(h)
   
    if density: 
        h.label="Density"
        hist.plot1d(h[dataset],density=True)
    else:
        hist.plot1d(h[dataset],density=False)

    if logy is not None:
        axs.set_yscale("log")
        axs.set_ylim(logy['min'],logy['max'])
    if ylim is not None:
        axs.set_ylim(ylim['min'],ylim['max'])
    if xlim is not None:
        axs.set_xlim(xlim['min'],xlim['max'])
    
    hep.cms.label(ax=axs)
    fig.savefig(outf)
    if close:
        plt.close(fig)
    return fig,axs

 
if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--loadhist', dest='loadhist', action='store_true',default = False, help='Load result from pickle')
    (options, args) = parser.parse_args()


    plots = [
        "gLLP_e","gLLP_pt","gLLP_eta","glepdPhi","gWPt"
    ]
    channels = ["muon","ele"]
    outdir = "./notebook/figures/"

    for channel in channels:
        
        if channel =="muon":
            #out = loadPickle("HNL_histograms_Mar15_muons_signal.pickle")
            out = loadPickle("HNL_histograms_Nov13_muons_signal.pickle")
            signals = [
             'HNL_muonType_mHNL1p0_pl10000_rwctau20000',
             'HNL_muonType_mHNL2p0_pl10000_rwctau5000',
             'HNL_muonType_mHNL4p0_pl1000',
             #'HNL_muonType_mHNL4p0_pl100_rwctau80',
            ]
        else:
            #out = loadPickle("HNL_histograms_Mar15_ele_signal.pickle")
            out = loadPickle("HNL_histograms_Nov13_ele_signal.pickle")
            signals = ['HNL_electronType_mHNL1p0_pl10000',
             'HNL_electronType_mHNL2p0_pl10000_rwctau5000',
             'HNL_electronType_mHNL4p0_pl1000',
             #'HNL_electronType_mHNL4p0_pl100_rwctau80'
        ]
    
        signals_bkg=signals+[     'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8']
        bkg = [     'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8']

        ## GEN plots
        regions = ["gLLP_csc","gLLP_dt"]
        for region in regions: 
            #plot(out['gLLP_e'].integrate("region",region), signals, outdir+"gLLP_e_%s_%s.pdf"%(channel,region))
            #plot(out['gLLP_pt'].integrate("region",region), signals, outdir+"gLLP_pt_%s_%s.pdf"%(channel,region))
            #plot(out['gLLP_eta'].integrate("region",region), signals, outdir+"gLLP_eta_%s_%s.pdf"%(channel,region),True,None,{"min":None,"max":1.2})
            #plot(out['glepdPhi'].integrate("region",region), signals, outdir+"glepdPhi_%s_%s.pdf"%(channel,region))
            pass

        #plot(out['gWPt'], signals_bkg, outdir+"gWPt_%s.pdf"%(channel))


        ## Reco plots
        regions = ["noselection","gLLP_csc","gLLP_dt"]
        for region in regions: 
            #plot(out['metXYCorr'].integrate("region",region), signals_bkg, outdir+"MET_%s_%s.pdf"%(channel,region),True,{"min":1e-5,"max":None})
            pass
        
        #plot(out['nLeptons'], signals_bkg, outdir+"nLeptons_%s.pdf"%(channel),True,{"min":1e-6,"max":10})
        #plot(out['nJets'], signals_bkg, outdir+"nJets_%s.pdf"%(channel))
        #plot(out['jetMet_dPhi'], signals_bkg, outdir+"jetMet_dPhi_%s.pdf"%(channel))

        if channel == "muon":
            #plot(out['muPt'], signals_bkg, outdir+"muPt_%s.pdf"%(channel))
            fig,ax = plot(out['muEta'], signals_bkg, outdir+"muEta_%s.pdf"%(channel))
            ax.set_ylim(0,0.7)
            fig.savefig(outdir+"muEta_muon.pdf")
        else:
            #plot(out['elePt'], signals_bkg, outdir+"elePt_%s.pdf"%(channel))
            fig,ax = plot(out['eleEta'], signals_bkg, outdir+"eleEta_%s.pdf"%(channel))
            ax.set_ylim(0,0.7)
            fig.savefig(outdir+"eleEta_ele.pdf")
            
        # cluster eff
        for signal in signals:
            ## CSC
            hr = out['llp_cls_eff_r']
            hz = out['llp_cls_eff_z']
            fname = outdir+"cls_eff_CSCrz_%s_%s.pdf"%(channel,rename(signal))
            plotEffrz(hr,hz,[signal],fname,"CSC")
            ## CSC
            hr = out['llp_cls_dt_eff_r']
            hz = out['llp_cls_dt_eff_z']
            fname = outdir+"cls_eff_DTrz_%s_%s.pdf"%(channel,rename(signal))
            plotEffrz(hr,hz,[signal],fname,"DT")
        
        #plotEff_e(out["llp_cls_eff_e"],signals, outdir+"cls_eff_E_csc_%s.pdf"%channel)
        #plotEff_e(out["llp_cls_dt_eff_e"],signals, outdir+"cls_eff_E_dt_%s.pdf"%channel)

            
        ## ABCD plots
        region = "ABCD"
        h = out['dphi_cluster_csc'].integrate("region",region)
        #plot(h.project("ClusterSize","dataset").rebin("ClusterSize",5), signals, outdir+"clsSize_csc_%s_%s.pdf"%(channel,region),True,{"min":1e-6,"max":10})
        #plot(h.project("dphi_lep","dataset").rebin("dphi_lep",2)   , signals, outdir+"dphi_lep_csc_%s_%s.pdf"%(channel,region))
        #plot(h.project("dphi_MET","dataset").rebin("dphi_MET",2)   , signals, outdir+"dphi_MET_csc_%s_%s.pdf"%(channel,region))

        region = "ABCD_dt"
        h = out['dphi_cluster_dt'].integrate("region",region)
        #plot(h.project("ClusterSize","dataset").rebin("ClusterSize",5), signals, outdir+"clsSize_dt_%s_%s.pdf"%(channel,region),True,{"min":1e-6,"max":10})
        #plot(h.project("dphi_lep","dataset").rebin("dphi_lep",2)      , signals, outdir+"dphi_lep_dt_%s_%s.pdf"%(channel,region))
        #plot(h.project("dphi_MET","dataset").rebin("dphi_MET",2)      , signals, outdir+"dphi_MET_dt_%s_%s.pdf"%(channel,region))



