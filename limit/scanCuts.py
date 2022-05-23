
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

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('-f', dest='fin', default = "", help='pickle name ')
    parser.add_option('--muon', dest='muon', action='store_true',default = False, help='make muon datacard')
    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='print cmd only')
    parser.add_option('--writeYields', dest='writeYields', action='store_true',default = False, help='write yield fields')
    parser.add_option('--combine', dest='combine', action='store_true',default = False, help='combine CSC and DT')
    (options, args) = parser.parse_args()

    #outdir = "../combine/dt_datacards/"
    outdir = "./combine/HNL_datacards/"
    suffix = ""

    bkg_unc = {
#        "bkg_sys_A":[0.5,0,0 ,0  ],
        "bkg_sys_D":[0  ,0,0 ,0.5],
    }

    sig_unc = {
        "clus":{
            "ggH":[0.15 ,0.15, 0.15 ,0.15],
        },
        "lumi":{
            "ggH":[0.018 ,0.018, 0.018 ,0.018],
        },
    }

    dphi_lepAxis = np.linspace(0,np.pi,31)[1:-2]
    if options.muon:
        cut = {"CSC":(220,dphi_lepAxis[-10],0.7), "DT":(110,dphi_lepAxis[-8],0.7)}
    else:
        cut = {"CSC":(160,dphi_lepAxis[-2],0.7), "DT":(100,dphi_lepAxis[-10],0.7)}

    #outdir = "./combine/HNL_datacards/scan/"
    #outdir = "./combine/HNL_datacards/scan_v2_fixed/"
    outdir = "./combine/HNL_datacards/scan_v3/" ## toysFreq, 50% 
    scanPoints = 0
    for channel in ["ele","muon"]:
        if channel == "muon": options.muon =True
        else: options.muon =False

        for region in ["CSC","DT"]:
            if not (channel=="muon" and region=="DT"): continue
            if channel=="ele" and region =="DT":
                sizeCuts = np.array([80,90,100,110,120,130])
            elif channel=="muon" and region =="DT":
                #sizeCuts = np.array([80,90,100,110,120,130])
                sizeCuts = np.array([140])
            elif channel=="ele" and region =="CSC":
                sizeCuts = np.array([140,160,180,200,220])
            elif channel=="muon" and region =="CSC":
                sizeCuts = np.array([160,180,200,220,240])
                #sizeCuts = np.array([230,240,250])

            #dphi_lepCuts = dphi_lepAxis[(dphi_lepAxis>1.0)& (dphi_lepAxis<2.5)]
            #dphi_METCuts = np.array([0.6,0.7,0.8,0.9,1.0])
            #dphi_lepCuts = dphi_lepAxis[(dphi_lepAxis>1.0)& (dphi_lepAxis<3.0)][::3]
            dphi_lepCuts = np.array([1.04719755, 1.36135682, 1.67551608, 1.98967535, 2.30383461,  2.61799388, 2.72271363, 2.82743339, 2.93215314])
            #dphi_METCuts = np.array([1.0,1.1,1.2,1.3,1.4])
            dphi_METCuts = np.array([0.7,1.0,None])
            #dphi_METCuts = np.array([1.0])
            #dphi_METCuts = np.array([1.0,1.4,1.6,2.0])
            for cut in cartesian_product(sizeCuts,dphi_lepCuts,dphi_METCuts):
               
                scanPoints+=1
                cutname = "s%s_lep%s_met%s"%(int(cut[0]),np.round(cut[1],2),str(cut[2]))
                outf = outdir + "yields_%s_%s_%s.json"%(channel,region,cutname)

                print(cut,cutname) 
                if region =="CSC":
                    cut = {"CSC":(cut[0],cut[1],cut[2]), "DT":(100,dphi_lepAxis[-10],dphi_lepAxis[-1])}
                else:
                    cut = {"CSC":(160,dphi_lepAxis[-2],dphi_lepAxis[-1]), "DT":(cut[0],cut[1],cut[2])}

                if options.writeYields:
                    from writeABCD import writeYields
                    writeYields(cut,options.muon,outf,False)
                else: 
                    with open(outf,'r') as f:
                        data = json.load(f)
                    bkg_rate_CSC = data['bkg']["CSC"]
                    bkg_rate_DT = data['bkg']["DT"]
                   
                    for name,signal in data.items():
                        if not ("HNL" and "mHNL4p0_pl1000") in name: continue
                        suffix = cutname
                        name = name +"_"+ suffix
                        #norm = 1
                        norm = signal["norm"]
                        from writeABCD import make_datacard_2sig 
                        if region == "CSC":    
                            sigRate = {"ggH":np.array(signal["CSC"])/norm }
                            obs = bkg_rate_CSC
                            obs[-1] = bkg_rate_CSC[0]*bkg_rate_CSC[2]/bkg_rate_CSC[1]  ## force D=A*C/B
                            make_datacard_2sig(outdir,name+"_CSC", sigRate, norm, bkg_rate_CSC, obs, bkg_unc,  sig_unc)
    
                            csc_limit = "combine -M AsymptoticLimits {odir}{name}_CSC.txt -n _{name}_CSC --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
                            print(csc_limit)
                            if not options.dryRun:
                                os.system(csc_limit)
                                os.system("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_CSC",outdir))
                        if region == "DT":    
                            dt_limit = "combine -M AsymptoticLimits {odir}{name}_DT.txt -n _{name}_DT --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
                            print(dt_limit)
                            sigRate = {"ggH":np.array(signal["DT"]) /norm}
                            obs = bkg_rate_DT
                            obs[-1] = bkg_rate_DT[0]*bkg_rate_DT[2]/bkg_rate_DT[1]
                            make_datacard_2sig(outdir,name+"_DT", sigRate, norm, bkg_rate_DT, obs, bkg_unc,  sig_unc)
    
                            if not options.dryRun:
                                os.system(dt_limit)
                                os.system("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_DT",outdir))
            
                        if options.combine:
                            cmd = "python combination.py CSC={odir}{name}_CSC.txt DT={odir}{name}_DT.txt {odir}{name}_comb.txt".format(name=name,odir=outdir)
                            print(cmd)
                            if not options.dryRun:
                                os.system(cmd)
                            cmd = "combine -M AsymptoticLimits {odir}{name}_comb.txt -n _{name}_comb --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
                            print(cmd)
                            if not options.dryRun:
                                os.system(cmd)
                                os.system("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_comb",outdir))
            print("Total scan points = ",scanPoints)
