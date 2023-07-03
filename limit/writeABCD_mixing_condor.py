import numpy as np
import os
import sys
sys.path.insert(0,"../")
import pickle
import json
from optparse import OptionParser

from writeABCD import make_datacard_2sig 

def Run(cmd,test=False):
    print(cmd)
    if not test: os.system(cmd)
    return
def RunLimit(outdir,name,tag1,tag2=None,tag3=None,test=False):
    if tag2 is None:
        #cmd = "combine -M AsymptoticLimits {odir}{name}_{tag1}.txt -n _{name}_{tag1} --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,tag1=tag1,odir=outdir,norm=1)
        cmd = "combine -M AsymptoticLimits {odir}{name}_{tag1}.txt -n _{name}_{tag1} --setParameters norm={norm} --freezeParameter norm ".format(name=name,tag1=tag1,odir=outdir,norm=1)
        Run(cmd,test)
        Run("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_"+tag1,outdir),test)
    else:
        cmd = "python combination.py {tag1}={odir}{name}_{tag1}.txt {tag2}={odir}{name}_{tag2}.txt {odir}{name}_{tag3}.txt".format(name=name,tag1=tag1,tag2=tag2,tag3=tag3,odir=outdir)
        Run(cmd,test)
        #cmd = "combine -M AsymptoticLimits {odir}{name}_{tag3}.txt -n _{name}_{tag3} --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,tag3=tag3,odir=outdir,norm=1)
        cmd = "combine -M AsymptoticLimits {odir}{name}_{tag3}.txt -n _{name}_{tag3} --setParameters norm={norm} --freezeParameter norm ".format(name=name,tag3=tag3,odir=outdir,norm=1)
        Run(cmd,test)
        Run("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_"+tag3,outdir),test)
    return

def makeAllcards_mixing(f_yield,muon,outdir="./combine/HNL_datacards/",test=False,fe=1.0,fmu=0.0,ftau=0.0,mass=1.0):
    mixed_name = "HNL_mixed-fe{}-fmu{}-ftau{}_".format(fe,fmu,ftau).replace(".","p")
    print(mixed_name)
    from collections import OrderedDict
    if not os.path.exists(outdir):
        print("mk dir = ", outdir)
        os.mkdir(outdir)
    norm = 100.0    # 1% BR

    with open(f_yield,'r') as f:
        data = json.load(f)
        if "DT_MB2" in data['bkg'].keys():
            hasMB2=True
        else:
            hasMB2 = False

    if muon:
        bkg_proc_CSC= {    "Zmumu_CSC":[-1,-1,-1,4.1]} ## 4.8% TF, 200 cut, full lumi}
        #bkg_unc_CSC = {    "Zmumu_CSC":[-1,-1,-1,0.26]    }
        bkg_unc_CSC = {    "Zmumu_CSC":[-1,-1,-1,0.30]    } #new systematics
        ## single-bin DT with total Zmumu
        bkg_proc_DT     = {           "Zmumu_DT"        :[-1,-1,-1,7.2]    }
        bkg_unc_DT         = {           "Zmumu_DT"        :[-1,-1,-1,0.42]    }
        ## two-bin DT with Zmumu
        #bkg_proc_DT_MB2 = {           "Zmumu_DT_MB2"        :[-1,-1,-1,7]      }# 36% TF @150         
        bkg_proc_DT_MB2 = {           "Zmumu_DT_MB2"        :[-1,-1,-1,8.2]      }# 36% TF @150,full lumi         
        #bkg_unc_DT_MB2  = {           "Zmumu_DT_MB2_sys_D"  :[0 ,0 ,0 ,0.42]   }    # 0.42=3/7         
        #bkg_unc_DT_MB2  = {           "Zmumu_DT_MB2_sys_D"  :[0 ,0 ,0 ,0.95]   }    # 0.95 = 6.9/7.2
        bkg_unc_DT_MB2  = {           "Zmumu_DT_MB2_sys_D"  :[0 ,0 ,0 ,0.90]   }    # 0.95 = 7.41/8.22
        #bkg_proc_DT_MB34= {           "Zmumu_DT_MB34"       :[-1,-1,-1,0.2,0.2]}
        #bkg_proc_DT_MB34= {           "Zmumu_DT_MB34"       :[-1,-1,-1,0.05,0.05]} # constant TF
        bkg_proc_DT_MB34= {           "Zmumu_DT_MB34"       :[-1,-1,-1,0.06,0.06]} # constant TF, full lumi
        bkg_unc_DT_MB34 = {           "Zmumu_DT_MB34_sys_D" :[0 ,0 ,0 ,1]      }# 100%  unc         
    else:
        bkg_proc_CSC= {}
        bkg_proc_DT= {}
        bkg_unc = {}
        bkg_unc_CSC = {}
        bkg_unc_DT = {}


    sig_unc_DT = OrderedDict({
        "lumi"            :{       "HNL":[0.016,0.016, 0.016 ,0.016],    },
        "pileup"          :{       "HNL":[0.01,0.01, 0.01 ,0.01],    },
        "Wxsec"           :{       "HNL":[0.038,0.038, 0.038 ,0.038],        },
        "Wpt"             :{       "HNL":[0.016,0.016, 0.016 ,0.016],        },
        "JES"             :{       "HNL":[0.02,0.02, 0.02 ,0.02],        },     ## TODO update
        "dt_clusterSyst"  :{       "HNL":[0.16,0.16, 0.16 ,0.16],    },
        #"dt_rpcSyst"      :{       "HNL":[0.053,0.053, 0.053 ,0.053],    },
        "dt_MB1Syst"      :{       "HNL":[0.074,0.074, 0.074 ,0.074],    },
        "dt_mixingSyst"  :{       "HNL":[0.10,0.10, 0.10 ,0.10],    },
    })
    sig_unc_CSC = OrderedDict({
        "lumi"            :{       "HNL":[0.016,0.016, 0.016 ,0.016],    },
        "pileup"          :{       "HNL":[0.01,0.01, 0.01 ,0.01],    },         ## TODO update
        "Wxsec"           :{       "HNL":[0.038,0.038, 0.038 ,0.038],        },
        "Wpt"             :{       "HNL":[0.016,0.016, 0.016 ,0.016],        },
        "JES"             :{       "HNL":[0.02,0.02, 0.02 ,0.02],        },     ## TODO update
        "csc_clusterSyst" :{       "HNL":[0.13 ,0.13 , 0.13  ,0.13 ],        },
        "csc_muonVeto"    :{       "HNL":[0.045 ,0.045, 0.045 ,0.045],        },
        "csc_jetVeto"     :{       "HNL":[0.0006 ,0.0006 , 0.0006  ,0.0006 ],        },
        "csc_rechitVeto"  :{       "HNL":[0.001 ,0.001 , 0.001  ,0.001 ],        },
        "csc_cut_based_id":{       "HNL":[0.051 ,0.051 , 0.051  ,0.051 ],        },
        "csc_time"        :{       "HNL":[0.009 ,0.009 , 0.009  ,0.009 ],        },
        "csc_time_spread" :{       "HNL":[0.028 ,0.028 , 0.028  ,0.028 ],        },
        "csc_readout"     :{       "HNL":[0.01 ,0.01 , 0.01  ,0.01 ],        },
        "csc_mixingSyst"  :{       "HNL":[0.10,0.10, 0.10 ,0.10],    },
    })
    if muon:
        mu_sys = { 
            "muTrig": {"HNL":[0.01,0.01, 0.01 ,0.01]},
            "muID"  : {"HNL":[0.003,0.003, 0.003 ,0.003]},
            "muISO" : {"HNL":[0.006,0.006, 0.006 ,0.006]},
        }
        sig_unc_CSC.update( mu_sys)
        sig_unc_DT.update( mu_sys)
    else:
        ele_sys = { 
            "eleTrig": {"HNL":[0.01,0.01, 0.01 ,0.01]},
            "eleReco" :{"HNL":[0.01,0.01, 0.01 ,0.01]},
        }
        sig_unc_CSC.update( ele_sys)
        sig_unc_DT.update( ele_sys)

    with open(f_yield,'r') as f:
        data = json.load(f)
        bkg_rate_CSC = {}       # {bkg: [a,b,c,d]}
        bkg_rate_CSC.update({"bkg":data['bkg']["CSC"]})          # {bkg: [a,b,c,d]}
        bkg_rate_CSC.update(bkg_proc_CSC)                        # add other bkg 
        bkg_rate_DT = {}       # {bkg: [a,b,c,d]}
        bkg_rate_DT.update({"bkg":data['bkg']["DT"]})           # {bkg: [a,b,c,d]}
        bkg_rate_DT.update(bkg_proc_DT)                         # add other bkg 
        if hasMB2 & muon:
            bkg_rate_DT_MB2 = {} 
            bkg_rate_DT_MB2.update({"bkg":data['bkg']["DT_MB2"]})
            bkg_rate_DT_MB2.update(bkg_proc_DT_MB2)
            bkg_rate_DT_MB34 = {} 
            bkg_rate_DT_MB34.update({"bkg":data['bkg']["DT_MB34"]})
            bkg_rate_DT_MB34.update(bkg_proc_DT_MB34)
            
        for name,signal in data.items():
            if "bkg" in name: continue
            if mixed_name not in name:  continue
            print("Found yields:", name)
            norm = 1
            #norm = signal["norm"]

            sigRate = {"HNL":np.array(signal["CSC"])/norm }
            obs = bkg_rate_CSC['bkg']
            #obs[-1] = bkg_rate_CSC['bkg'][0]*bkg_rate_CSC['bkg'][2]/bkg_rate_CSC['bkg'][1]  ## force D=A*C/B
            make_datacard_2sig(outdir,name+"_CSC", sigRate, norm, bkg_rate_CSC, obs, bkg_unc_CSC,  sig_unc_CSC)
            sigRate = {"HNL":np.array(signal["DT"]) /norm}
            obs = bkg_rate_DT['bkg']
            #obs[-1] = bkg_rate_DT['bkg'][0]*bkg_rate_DT['bkg'][2]/bkg_rate_DT['bkg'][1]  ## force D=A*C/B
            make_datacard_2sig(outdir,name+"_DT", sigRate, norm, bkg_rate_DT, obs, bkg_unc_DT,  sig_unc_DT)
            if hasMB2 & muon:
                obs = bkg_rate_DT_MB2['bkg']
                sigRate = {"HNL":np.array(signal["DT_MB2"]) /norm}
                #obs[-1] = bkg_rate_DT_MB2['bkg'][0]*bkg_rate_DT_MB2['bkg'][2]/bkg_rate_DT_MB2['bkg'][1]  ## force D=A*C/B
                make_datacard_2sig(outdir,name+"_DT_MB2", sigRate, norm, bkg_rate_DT_MB2, obs, bkg_unc_DT_MB2,  sig_unc_DT)
                obs = bkg_rate_DT_MB34['bkg']
                sigRate = {"HNL":np.array(signal["DT_MB34"]) /norm}
                #obs[-1] = bkg_rate_DT_MB34['bkg'][0]*bkg_rate_DT_MB34['bkg'][2]/bkg_rate_DT_MB34['bkg'][1]  ## force D=A*C/B
                make_datacard_2sig(outdir,name+"_DT_MB34", sigRate, norm, bkg_rate_DT_MB34, obs, bkg_unc_DT_MB34,  sig_unc_DT)
   
                    
            RunLimit(outdir,name,"CSC",None,None,test)
            if hasMB2:
                #RunLimit(outdir,name,"DT_MB2",None,None,test)
                #RunLimit(outdir,name,"DT_MB34",None,None,test)
                RunLimit(outdir,name,"DT_MB2","DT_MB34","DTcomb",test)
                pass
            else:
                RunLimit(outdir,name,"DT",None,None,test)

            #Combination
            if hasMB2 & muon:
                RunLimit(outdir,name,"CSC","DTcomb","comb",test)       # CSC+MB2+MB34
                #RunLimit(outdir,name,"CSC","DT_MB34","combNoMB2",test) # CSC+MB34
                #RunLimit(outdir,name,"CSC","DT","combSingleDTBin",test)            # CSC+DT
                #RunLimit(outdir,name,"CSC","DT_MB34","comb",test) # CSC+MB34
                pass
            else:
                RunLimit(outdir,name,"CSC","DT","comb",test)
            
    return 

def combine_channels(f_yield,outdir_ele,outdir_mu,outdir_comb,fe=1.0,fmu=0.0,ftau=0.0,test=False):
    mixed_name = "HNL_mixed-fe{}-fmu{}-ftau{}_".format(fe,fmu,ftau).replace(".","p")
    print(mixed_name)
    with open(f_yield,'r') as f:
        data = json.load(f)
    if not os.path.exists(outdir_comb):
        print("mk dir = ", outdir_comb)
        os.mkdir(outdir_comb)

    for name,signal in data.items():
            if "bkg" in name: continue
            if mixed_name not in name:  continue
            if os.path.exists("{odir_ele}{name}_comb.txt".format(odir_ele=outdir_ele,name=name)) and os.path.exists("{odir_mu}{name}_comb.txt".format(odir_mu=outdir_mu,name=name)):
                cmd="python combination.py Ele={odir_ele}{name}_comb.txt Muon={odir_mu}{name}_comb.txt {odir_comb}{name}_comb.txt".format(odir_ele=outdir_ele,name=name,odir_mu=outdir_mu,odir_comb=outdir_comb)
                Run(cmd,test)

                cmd = "combine -M AsymptoticLimits {odir}{name}_comb.txt -n _{name}_comb --setParameters norm={norm} --freezeParameter norm ".format(name=name,odir=outdir_comb,norm=1)
                            
                Run(cmd,test)

                Run("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_comb",outdir_comb),test)
if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option('-f', dest='fin', default = "", help='yield json name ')

    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='dryRun')
    parser.add_option('--eleJson', dest='eleJson', default = "", help='json files for ele channels')
    parser.add_option('--muJson', dest='muJson', default = "", help='json files for muon channels')
    parser.add_option('--fe', dest='fe', default = "", help='fe mixing parameter')
    parser.add_option('--fmu', dest='fmu', default = "", help='fmu mixing parameter')
    parser.add_option('--ftau', dest='ftau', default = "", help='ftau mixing parameter')
    parser.add_option('--mass', dest='mass', default = "", help='HNL mass')


    (options, args) = parser.parse_args()
    outdir_ele = "./HNL_datacards_ele/"
    outdir_muon = "./HNL_datacards_muon/"
    outdir_comb = "./HNL_datacards_MuE/"

    
    if float(options.fe)>=0: 
        makeAllcards_mixing(options.eleJson,False,outdir_ele,options.dryRun,options.fe,options.fmu,options.ftau,options.mass) 
    if float(options.fmu)>=0: 
        makeAllcards_mixing(options.muJson,True,outdir_muon,options.dryRun,options.fe,options.fmu,options.ftau,options.mass)
    if float(options.ftau)>=0:
        combine_channels(options.muJson,outdir_ele,outdir_muon,outdir_comb,options.fe,options.fmu,options.ftau)
        ## Write functions to combine ele/muon output
       
        pass
