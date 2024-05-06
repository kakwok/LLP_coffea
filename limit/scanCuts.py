
import os
import sys
import json
import numpy as np
from optparse import OptionParser
import awkward as ak
import json

def cartesian_product(*arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)

def getBkgUnc(cut,options):
    if options.muon:
        # must use convention :     procName_uncName
        bkg_unc = {
            #"bkg_sys_D":[0  ,0 ,0 ,0.5],
            #"Zmumu_CSC_sys_D":[0 ,0 ,0 ,0.32],      
            #"Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.2 ],
            "Zmumu_CSC_sys_D":[0 ,0 ,0 ,0.29], ## 5.69% TF, 200 cut 
            #"Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.24], # 1.7% TF, 130, loose ID
            "Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.3], # 1.47% TF, 130, no loose ID
        }
        if cut['DT'][0]==100:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,3.70]}             
        if cut['DT'][0]==110:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,2.97]}             
        if cut['DT'][0]==120:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,2.76]}             
        if cut['DT'][0]==130:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,2.46]}             
        if cut['DT'][0]==140:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,2.03]}             
        if cut['DT'][0]==150:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,1.81]}             
        if cut['DT'][0]==160:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,1.72]}             
        if cut['DT'][0]==180:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,1.28]}             
        if cut['DT'][0]==200:  bkg_proc_DT= {"Zmumu_DT":[-1,-1,-1,0.98]}             
        if cut['CSC'][0]==100: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,26.32]    }
        if cut['CSC'][0]==140: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,11.36]    }
        if cut['CSC'][0]==160: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,7.97]     }
        if cut['CSC'][0]==180: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,5.54]     }
        if cut['CSC'][0]==200: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,4.25]     }
        if cut['CSC'][0]==220: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,2.91]     }
        if cut['CSC'][0]==240: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,2.20]     }
        if cut['CSC'][0]==260: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,1.72]     }
        if cut['CSC'][0]==280: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,1.67]     }
        if cut['CSC'][0]==300: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,1.50]     }
        if cut['CSC'][0]==320: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,1.00]     }
        if cut['CSC'][0]==340: bkg_proc_CSC= {  "Zmumu_CSC":[-1,-1,-1,0.84]     }

    else:
        bkg_proc_CSC= {}
        bkg_proc_DT= {}
        bkg_unc = {}
    return bkg_proc_CSC,bkg_proc_DT,bkg_unc

def getSigUnc(options):
    from collections import OrderedDict

    sig_unc_DT = OrderedDict({
        "lumi"            :{       "HNL":[0.016,0.016, 0.016 ,0.016],    },
        "pileup"          :{       "HNL":[0.01,0.01, 0.01 ,0.01],    },
        "Wpt"             :{       "HNL":[0.10,0.10, 0.10 ,0.10],        },
        "JES"             :{       "HNL":[0.02,0.02, 0.02 ,0.02],        },     ## TODO update
        "dt_clusterSyst"  :{       "HNL":[0.15,0.15, 0.15 ,0.15],    },
        "dt_rpcSyst"      :{       "HNL":[0.053,0.053, 0.053 ,0.053],    },
        "dt_MB1Syst"      :{       "HNL":[0.074,0.074, 0.074 ,0.074],    },
    })
    sig_unc_CSC = OrderedDict({
        "lumi"            :{       "HNL":[0.016,0.016, 0.016 ,0.016],    },
        "pileup"          :{       "HNL":[0.01,0.01, 0.01 ,0.01],    },         ## TODO update
        "Wpt"             :{       "HNL":[0.10,0.10, 0.10 ,0.10],        },     ## TODO update
        "JES"             :{       "HNL":[0.02,0.02, 0.02 ,0.02],        },     ## TODO update
        "csc_muonVeto"    :{       "HNL":[0.045 ,0.045, 0.045 ,0.045],        },
        "csc_jetVeto"     :{       "HNL":[0.00068 ,0.00068 , 0.00068  ,0.00068 ],        },
        "csc_rechitVeto"  :{       "HNL":[0.001 ,0.001 , 0.001  ,0.001 ],        },
        "csc_cut_based_id":{       "HNL":[0.051 ,0.051 , 0.051  ,0.051 ],        },
        "csc_time"        :{       "HNL":[0.009 ,0.009 , 0.009  ,0.009 ],        },
        "csc_time_spread" :{       "HNL":[0.028 ,0.028 , 0.028  ,0.028 ],        },
        "csc_clusterSyst" :{       "HNL":[0.035 ,0.035 , 0.035  ,0.035 ],        },
        "csc_readout"     :{       "HNL":[0.01 ,0.01 , 0.01  ,0.01 ],        },
    })
    if options.muon:
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
    return sig_unc_DT,sig_unc_CSC




if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('-f', dest='fin', default = "", help='pickle name ')
    parser.add_option('--muon', dest='muon', action='store_true',default = False, help='make muon datacard')
    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='print cmd only')
    parser.add_option('--writeYields', dest='writeYields', action='store_true',default = False, help='write yield fields')
    parser.add_option('--unblind', dest='unblind', action='store_true',default = False, help='do not blind bin D if true ')
    parser.add_option('--useOOT', dest='useOOT', action='store_true',default = False, help='useOOT for bkg ABCD prediction, use real ABC otherwise ')
    parser.add_option('--combine', dest='combine', action='store_true',default = False, help='combine CSC and DT')
    (options, args) = parser.parse_args()

    outdir = "./combine/HNL_datacards/"
    suffix = ""

    outdir = "./combine/HNL_datacards/scan_v4/" ## CSC,DT reopt, useOOT
    #outdir = "./combine/HNL_datacards/scan_v5/" ## CSC,DT reopt, use real ABC
    scanPoints = 0
    def scanNhit(channel,region):
        dphi_lepCuts = np.array([2.8])
        if  region =="DT":
            sizeCuts = np.array([100,110,120,130,140,150,160,180,200])
        elif  region =="CSC":
            sizeCuts = np.array([100,140,160,180,200,220,240,260,280,300,320,340])
        return dphi_lepCuts,sizeCuts

    def scanDphi(channel,region):
        if region=="CSC":       sizeCuts = np.array([200])
        if region=="DT":        sizeCuts = np.array([150])
        dphi_lepCuts = np.array([1.0,1.4,1.6,2.0,2.5,2.6,2.7,2.8,2.9,3.0])
        return dphi_lepCuts,sizeCuts


    print("output directory = ",outdir)
    for channel in ["ele","muon"]:
        if channel == "muon": options.muon =True
        else: options.muon =False

        dphi_METCuts = np.array([None])
        for region in ["CSC","DT"]:
            #if not (region=="DT"): continue
            #dphi_lepCuts,sizeCuts = scanNhit(channel,region)
            dphi_lepCuts,sizeCuts = scanDphi(channel,region)

            for cut in cartesian_product(sizeCuts,dphi_lepCuts,dphi_METCuts):
                
                scanPoints+=1
                cutname = "s%s_lep%s_met%s"%(int(cut[0]),np.round(cut[1],2),str(cut[2]))
                outf = outdir + "yields_%s_%s_%s.json"%(channel,region,cutname)

                print(cut,cutname) 
                if region =="CSC":
                    cut = {"CSC":(cut[0],cut[1],cut[2]), "DT":(130,cut[1],cut[2])}
                else:
                    cut = {"CSC":(160,cut[1],cut[2]), "DT":(cut[0],cut[1],cut[2])}

                if options.writeYields:
                    from writeABCD import writeYields
                    shifts={}
                    tauSignals=False
                    writeYields(cut,options.muon,outf,False,shifts,tauSignals,options.unblind,options.useOOT)
                else: 
                    with open(outf,'r') as f:
                        data = json.load(f)

                    bkg_proc_CSC,bkg_proc_DT,bkg_unc = getBkgUnc(cut,options)
                    sig_unc_DT,sig_unc_CSC           = getSigUnc(options)

                    bkg_rate_CSC = {}       # {bkg: [a,b,c,d]}
                    bkg_rate_CSC.update({"bkg":data['bkg']["CSC"]})          # {bkg: [a,b,c,d]}
                    bkg_rate_CSC.update(bkg_proc_CSC)                        # add other bkg 
                    bkg_rate_DT = {}       # {bkg: [a,b,c,d]}
                    bkg_rate_DT.update({"bkg":data['bkg']["DT"]})           # {bkg: [a,b,c,d]}
                    bkg_rate_DT.update(bkg_proc_DT)                         # add other bkg 
                   
                    for name,signal in data.items():
                        if not ("HNL" and "mHNL4p0_pl1000") in name: continue
                        suffix = cutname
                        name = name +"_"+ suffix
                        #norm = 1
                        norm = signal["norm"]
                        from writeABCD import make_datacard_2sig 
                        if region == "CSC":    
                            sigRate = {"HNL":np.array(signal["CSC"])/norm }
                            obs = bkg_rate_CSC['bkg']
                            obs[-1] = bkg_rate_CSC['bkg'][0]*bkg_rate_CSC['bkg'][2]/bkg_rate_CSC['bkg'][1]  ## force D=A*C/B
                            make_datacard_2sig(outdir,name+"_CSC", sigRate, norm, bkg_rate_CSC, obs, bkg_unc,  sig_unc_CSC)
    
                            csc_limit = "combine -M AsymptoticLimits {odir}{name}_CSC.txt -n _{name}_CSC --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
                            #csc_limit = "combine -M Significance {odir}{name}_CSC.txt -n _{name}_CSC --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq --expectSignal=1".format(name=name,odir=outdir,norm=1)
                            print(csc_limit)
                            if not options.dryRun:
                                os.system(csc_limit)
                                os.system("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_CSC",outdir))
                                #os.system("mv higgsCombine_%s.Significance.mH120.root %s"%(name+"_CSC",outdir))
                        if region == "DT":    
                            dt_limit = "combine -M AsymptoticLimits {odir}{name}_DT.txt -n _{name}_DT --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
                            #dt_limit = "combine -M Significance {odir}{name}_DT.txt -n _{name}_DT --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq --expectSignal=1".format(name=name,odir=outdir,norm=1)
                            print(dt_limit)
                            sigRate = {"HNL":np.array(signal["DT"]) /norm}
                            obs = bkg_rate_DT['bkg']
                            obs[-1] = bkg_rate_DT['bkg'][0]*bkg_rate_DT['bkg'][2]/bkg_rate_DT['bkg'][1]  ## force D=A*C/B
                            make_datacard_2sig(outdir,name+"_DT", sigRate, norm, bkg_rate_DT, obs, bkg_unc,  sig_unc_DT)
    
                            if not options.dryRun:
                                os.system(dt_limit)
                                os.system("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_DT",outdir))
                                #os.system("mv higgsCombine_%s.Significance.mH120.root %s"%(name+"_DT",outdir))
            
                        if options.combine:
                            cmd = "python combination.py CSC={odir}{name}_CSC.txt DT={odir}{name}_DT.txt {odir}{name}_comb.txt".format(name=name,odir=outdir)
                            print(cmd)
                            #if not options.dryRun:
                            #    os.system(cmd)
                            cmd = "combine -M AsymptoticLimits {odir}{name}_comb.txt -n _{name}_comb --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
                            print(cmd)
                            if not options.dryRun:
                                os.system(cmd)
                                os.system("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_comb",outdir))
        print("Total scan points = ",scanPoints)
