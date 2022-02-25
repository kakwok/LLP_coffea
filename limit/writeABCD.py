import numpy as np
import os
import sys
sys.path.insert(0,"../")
import pickle
from optparse import OptionParser

def make_datacard_2sig(outDataCardsDir,modelName,  signal_rate, norm, bkg_rate, observation, bkg_unc, bkg_unc_name, sig_unc):
    a,b,c,d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    #c1 = a/b
    #c1 = b/a        ## DT convention
    #c2 = c/b
    nSig = len(signal_rate.keys())
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# signal norm {0} \n'.format(norm))

    text_file.write('imax {0} \n'.format(4))
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')

    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \t chD \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n'.format(observation[0],observation[1],observation[2],observation[3]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig) +'\t chD '*(1+nSig) +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * 4 + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * 4 + '\n')
    rate_string = 'rate'
    for i in range(4):# 4 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')
    text_file.write('NA_val      extArg         {}    [const] \n'.format(a))
    text_file.write('NB_val      extArg          {}   [const] \n'.format(b))
    text_file.write('NC_val      extArg          {}   [const] \n'.format(c))
    text_file.write('NA_x        extArg          0       [-7,7] \n')
    text_file.write('NB_x        extArg          0       [-7,7] \n')
    text_file.write('NC_x        extArg          0       [-7,7] \n')
    text_file.write('NA  rateParam       chA     bkg      (@0*(1+1/sqrt(@0))**@1)       NA_val,NA_x \n')
    text_file.write('NB  rateParam       chB     bkg      (@0*(1+1/sqrt(@0))**@1)       NB_val,NB_x \n')
    text_file.write('NC  rateParam       chC     bkg      (@0*(1+1/sqrt(@0))**@1)       NC_val,NC_x \n')
    text_file.write('ND  rateParam       chD     bkg      (@0*@2/@1)                     NA,NB,NC ## D=A*C/B\n')
    text_file.write('norm  rateParam       *     ggH      1 \n')


def calABCD(listABCD,eff_nhit,eff_dphi):

    D = listABCD[3]
    C = (D / eff_nhit ) - D
    sumABCD = (C+D)/(eff_dphi)
    sumAB   = sumABCD - C - D
    B = sumAB*(D/(C+D))
    A = sumAB*(C/(C+D))
    listABCD[2] = np.round(C,3)
    listABCD[1] = np.round(B,3)
    listABCD[0] = np.round(A,3)
    return

    return
def calSigal():
    data = {
    #"HNL_muonType_mHNL1p0_pl10" :{'CSC': np.array([1313.72597472, 2083.29270286, 2387.01682551, 1588.47727032]), 'DT': np.array([1245.3276547 ,    0.        ,  514.90494597,    0.        ]), 'norm': 100},
    #"HNL_muonType_mHNL1p0_pl100" :{'CSC': np.array([  958.36080102,  2541.21051879, 35942.66029522, 15694.03822141]), 'DT': np.array([ 2516.41081957,  2687.36740861, 24510.1580307 , 20040.76760296]), 'norm': 100},
    #"HNL_muonType_mHNL1p0_pl1000": {'CSC': np.array([  28.82915809,   72.20854441, 1063.38960708,  482.29479901]), 'DT': np.array([ 55.81080617,  84.18083673, 980.13814336, 767.72333672]), 'norm': 10},
    #"HNL_muonType_mHNL1p0_pl2000": {'CSC':np.array([  7.5398423 ,  19.03590541, 286.4139182 , 127.38951007]), 'DT': np.array([ 14.8699539 ,  22.64145349, 266.42989897, 207.79041764]), 'norm': 10},
    #"HNL_muonType_mHNL1p0_pl5000": {'CSC':np.array([ 1.23971623,  3.1465619 , 48.01300414, 21.06873556]), 'DT': np.array([ 2.47262373,  3.78655542, 44.83523666, 34.87615238]), 'norm': 1},
    #"HNL_muonType_mHNL2p0_pl10": {'CSC': np.array([12.94954136,  0.        ,  8.14461429,  0.        ]), 'DT': np.array([ 2.89209086,  0.        , 12.84317809,  0.        ]), 'norm': 1},
    #"HNL_muonType_mHNL2p0_pl100": {'CSC': np.array([ 46.08049672,  65.31234612, 763.99850184, 382.83887428]), 'DT': np.array([ 66.71386847,  59.28390516, 373.94752874, 309.02756679]), 'norm': 1},
    #"HNL_muonType_mHNL2p0_pl1000": {'CSC': np.array([ 2.62177929,  6.13920546, 68.80480406, 22.27152984]), 'DT': np.array([ 4.84003038,  4.08741815, 62.3940202 , 39.06735284]), 'norm': 1},
    # from ct1m ###"HNL_muonType_mHNL2p0_pl15000": {'CSC': np.array([ 1.23607486,  2.94705472, 33.51721181, 10.69697318]), 'DT': np.array([ 2.34226239,  2.01891048, 30.98168154, 19.26501103]), 'norm': 1},
    ##from ct1m ###"HNL_muonType_mHNL2p0_pl2000": {'CSC': np.array([ 0.7164319 ,  1.72425638, 19.74654625,  6.25635821]), 'DT': np.array([ 1.37572843,  1.19858567, 18.42712308, 11.41849004]), 'norm': 1},
"HNL_muonType_mHNL2p0_pl10000":{'CSC': np.array([0.01535685, 0.08224103, 1.12456859, 0.30455309]), 'DT': np.array([0.06101222, 0.04465391, 0.90300775, 0.63003223]), 'norm': 1},
"HNL_muonType_mHNL2p0_pl2000" :{'CSC': np.array([ 0.35527838,  1.80945444, 25.12737918,  6.99742389]), 'DT': np.array([ 1.35857254,  1.04874446, 19.70587762, 13.900502  ]), 'norm': 1},
"HNL_muonType_mHNL2p0_pl5000" :{'CSC': np.array([0.06024815, 0.31859128, 4.37329474, 1.19273078]), 'DT': np.array([0.23706951, 0.17584326, 3.49085573, 2.44230647]), 'norm': 1},
    #"HNL_muonType_mHNL2p0_pl800": {'CSC': np.array([  3.92207907,   9.06671745, 100.43177256,  32.844178  ]), 'DT': np.array([ 7.10264082,  5.91173037, 89.76806317, 56.49657938]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl100": {'CSC': np.array([1.00382177, 0.88081391, 7.03360537, 5.26295055]), 'DT': np.array([0.91579247, 0.87776941, 1.97958913, 1.12285429]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl1000": {'CSC': np.array([0.11367208, 0.23031757, 3.90969984, 1.4149658 ]), 'DT': np.array([0.16994022, 0.33335487, 3.38758482, 1.77729034]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl1500": {'CSC': np.array([0.05674916, 0.12227352, 2.08696206, 0.7353985 ]), 'DT': np.array([0.08887886, 0.18265597, 1.8813124 , 0.97029353]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl200": {'CSC': np.array([ 0.89609229,  1.07542741, 13.89661744,  6.37708174]), 'DT': np.array([0.84830022, 1.07858285, 7.20642776, 4.46332663]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl2000": {'CSC': np.array([0.03388892, 0.07544135, 1.28896805, 0.44779054]), 'DT': np.array([0.05441054, 0.11464907, 1.18446408, 0.60551109]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl500": {'CSC': np.array([0.32837435, 0.56626363, 9.24723856, 3.59179129]), 'DT': np.array([0.43360458, 0.75001888, 7.07675428, 3.89266032]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl800": {'CSC': np.array([0.16322816, 0.31673343, 5.34028832, 1.96939856]), 'DT': np.array([0.23616548, 0.44824362, 4.48827116, 2.38402303]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl50": {'CSC': np.array([0.52972903, 0.41909375, 1.16959573, 1.26086066]), 'DT': np.array([0.38047228, 0.29476375, 0.12242464, 0.05909414]), 'norm': 1},
    #"HNL_muonType_mHNL4p0_pl80": {'CSC': np.array([0.9071468 , 0.73450322, 4.61030034, 3.84690764]), 'DT': np.array([0.71360755, 0.67620385, 0.99318524, 0.5558966 ]), 'norm': 1},
    }
    #data={
    #    "HNL_muonType_mHNL1p0_pl10":{
    #                "CSC" :[1313.72597472, 2083.29270286, 2387.01682551, 1588.47727032],
    #                "DT"  :[1245.3276547 ,    0.         , 514.90494597 ,   0.        ],
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL1p0_pl100":{
    #                "CSC" :[  958.36080102 , 2541.21051879, 35942.66029522 ,15694.03822141],
    #                "DT"  :[ 2516.41081957 , 2687.36740861, 24510.1580307  ,20040.76760296],
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL2p0_pl100":{
    #                "CSC" :[ 46.08049672,  65.31234612, 763.99850184, 382.83887428],
    #                "DT"  :[ 66.71386847,  59.28390516, 373.94752874, 309.02756679],
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL2p0_pl1000":{
    #                "CSC" :[ 2.62177929 , 6.13920546 ,68.80480406, 22.27152984],
    #                "DT"  :[ 4.84003038 , 4.08741815 ,62.3940202 , 39.06735284],
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL2p0_pl1500":{
    #                "CSC" :[ 1.23607486 , 2.94705472 ,33.51721181, 10.69697318],
    #                "DT"  :[ 2.34226239 , 2.01891048 ,30.98168154, 19.26501103],
    #                "norm":1
    #    },

    #    "HNL_muonType_mHNL4p0_pl50":{
    #                "CSC" :[0.52972903 ,0.41909375, 1.16959573 ,1.26086066] ,
    #                "DT"  :[0.38047228 ,0.29476375, 0.12242464 ,0.05909414] ,
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL4p0_pl100":{
    #                "CSC" :[1.00382177 ,0.88081391 ,7.03360537 ,5.26295055] ,
    #                "DT"  :[0.91579247 ,0.87776941 ,1.97958913 ,1.12285429] ,
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL4p0_pl1000":{
    #                "CSC" :[0.11367208, 0.23806195, 4.20743326, 1.46618138] ,
    #                "DT"  :[0.16994022, 0.33335487, 3.38758482, 1.77729034] ,
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL7p0_pl100":{
    #                "CSC" :[0.03804715 ,0.02123373 ,0.05201338 ,0.03323407] ,
    #                "DT"  :[0.0135374  ,0.01267536 ,0.01290687 ,0.00190483] ,
    #                "norm":1
    #    },
    #    "HNL_muonType_mHNL7p0_pl1000":{
    #                "CSC" :[0.00872521, 0.01415912 ,0.26534225, 0.11173698],
    #                "DT"  :[0.01352161, 0.02273869 ,0.22597267, 0.12207775],
    #                "norm":1
    #    },

    #}
    #for name,signal in data.items():
    #    calABCD(signal["MB2"],signal["nhit"],signal["dphi"])
    #    calABCD(signal["MB34"],signal["nhit"],signal["dphi"])
    for d in data:
        print(d,data[d])
    return data

def toInt(h):
    return list(h.values(True).values())[0]

def predIntimeFromOOT(h,size,dphi_lep,dphi_met,isSignal=False,kfactor=0.25,lumi=137/23.):
    # kfactor=0.25
    # lumi = 137/23.
    
    cut1 = slice(size,None)
    cut2 = slice(dphi_lep,None)
    cut3 = slice(None,dphi_met)
    A = h.integrate("ClusterSize",slice(size,None)).integrate("dphi_lep",slice(None,dphi_lep)).integrate("dphi_MET",cut3)
    B = h.integrate("ClusterSize",slice(None,size)).integrate("dphi_lep",slice(None,dphi_lep)).integrate("dphi_MET",cut3)
    C = h.integrate("ClusterSize",slice(None,size)).integrate("dphi_lep",slice(dphi_lep,None)).integrate("dphi_MET",cut3)
    D = h.integrate("ClusterSize",slice(size,None)).integrate("dphi_lep",slice(dphi_lep,None)).integrate("dphi_MET",cut3)
    
    if not isSignal:
        N_evts    = np.array([toInt(A)[0],toInt(B)[0],toInt(C)[0],toInt(D)[0]]) * kfactor * lumi
        N_evts_unc = np.sqrt(np.array([toInt(A)[1],toInt(B)[1],toInt(C)[1],toInt(D)[1]])) * kfactor * lumi        
    else:
        N_evts     = np.array([toInt(A)[0],toInt(B)[0],toInt(C)[0],toInt(D)[0]]) 
        N_evts_unc = np.sqrt(np.array([toInt(A)[1],toInt(B)[1],toInt(C)[1],toInt(D)[1]]) )
 
    return N_evts, N_evts_unc

def loadhist():
    from coffea import hist
    
    #with open('../HNL_histograms_Feb18_muons_signal.pickle','rb') as f:                
    with open('../HNL_histograms_Feb23_muons_signal.pickle','rb') as f:                
        out = pickle.load(f)
   
    lumi = 137.0
    
    for k,h in out.items():
        if (type(h)!=hist.Hist): continue
        h.scale({ d: lumi for d in h.identifiers("dataset")} , axis="dataset") 

    signalNames = [ s.name for s in out['dphi_cluster_csc'].identifiers("dataset")]   
    data = {} 
    for signal_name in signalNames:
        h = out['dphi_cluster_csc']
        signal = h.integrate("dataset",signal_name).integrate("region","ABCD")
        hdt = out['dphi_cluster_dt']
        signal_dt = hdt.integrate("dataset",signal_name).integrate("region","ABCD_dt")

        dphi_lepCuts = np.linspace(0,np.pi,31)[1:-2]
        ### NOT USED FOR SIGNAL
        lumi = 1
        kfactor=1     ## Muon, CSC InT/OOT 

        sizeCut = 220
        dphi_lepCut = dphi_lepCuts[-10] ###.1.989675
        dphi_METCut = 0.7
        dphi_lepCuts = np.linspace(0,np.pi,31)[1:-2]

        #print(predIntimeFromOOT(signal,sizeCut,dphi_lepCut,dphi_METCut,True,kfactor,lumi))
        CSC,CSC_unc = predIntimeFromOOT(signal,sizeCut,dphi_lepCut,dphi_METCut,True,kfactor,lumi)

        kfactor=0.9

        sizeCut = 110
        dphi_lepCut = dphi_lepCuts[-8] ###. 2.199115
        dphi_METCut = 0.7
        ## s =5.3, 
        #print(predIntimeFromOOT(signal_dt,sizeCut,dphi_lepCut,dphi_METCut,True,kfactor,lumi))
        DT,DT_unc = predIntimeFromOOT(signal_dt,sizeCut,dphi_lepCut,dphi_METCut,True,kfactor,lumi)
        if "rwctau" in signal_name:
            ct = signal_name.split("_")[-1].replace("rwctau","pl")
            sample_name = ("_".join(signal_name.split("_")[:-2]+[ct]))
        else:
            sample_name = signal_name 
        if "2p0" in sample_name:
            print(signal_name," CSC = " ,CSC)
            print(signal_name," CSCunc = " ,CSC_unc)
            print(signal_name," DT = " ,DT)
            print(signal_name," DTunc = ", DT_unc)
        data[sample_name] = {"CSC":CSC,"DT":DT,"norm":1}
              
    return data 

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--loadhist', dest='loadhist', action='store_true',default = False, help='Load result from pickle')
    (options, args) = parser.parse_args()

    #outdir = "../combine/dt_datacards/"
    outdir = "./combine/HNL_datacards/"
    norm = 100.0    # 1% BR
    suffix = ""

    
    bkg_rate_CSC_mu = [2.97826087, 856.25, 2114.56522, 1.48913043]
    bkg_rate_DT_mu =[  16.0826087,  1088.25652174 ,1688.67391304 ,   5.36086957] 
    bkg_unc_CSC_mu  = [0.3]
    bkg_unc_DT_mu  = [0.3]
    bkg_unc_name = ["unc1"]
    sig_unc = {
        "clus":{
            "ggH":[0.15 ,0.15, 0.15 ,0.15],
        },
        "lumi":{
            "ggH":[0.018 ,0.018, 0.018 ,0.018],
        },
    }

    if options.loadhist:
        data = loadhist()
        for k,v in data.items():
            print(k,v)
    else:
        data = calSigal()
        for name,signal in data.items():
            if suffix:
                name = name+suffix
            #norm = 1
            norm = signal["norm"]
            sigRate = {"ggH":np.array(signal["CSC"])/norm }
            obs = bkg_rate_CSC_mu
            make_datacard_2sig(outdir,name+"_CSC", sigRate, norm, bkg_rate_CSC_mu, obs, bkg_unc_CSC_mu, bkg_unc_name, sig_unc)
            sigRate = {"ggH":np.array(signal["DT"]) /norm}
            obs = bkg_rate_DT_mu
            make_datacard_2sig(outdir,name+"_DT", sigRate, norm, bkg_rate_DT_mu, obs, bkg_unc_DT_mu, bkg_unc_name, sig_unc)
    
            #sigRate = {"ggH":np.array(signal["MB34"]) * 100 /norm }
            #sigRate = {"ggH":np.array(signal["MB34"])}
            #obs = bkg_rate_MB34
            #make_datacard_2sig(outdir,name+"_MB34", sigRate, norm, bkg_rate_MB34, obs, bkg_unc_MB34, bkg_unc_name, sig_unc)
    
            cmd = "python combination.py CSC={odir}{name}_CSC.txt DT={odir}{name}_DT.txt {odir}{name}_comb.txt".format(name=name,odir=outdir)
            print(cmd)
            os.system(cmd)
            cmd = "combine -M AsymptoticLimits {odir}{name}_comb.txt -n _{name}_comb --setParameters norm={norm} --freezeParameter norm -t -1".format(name=name,odir=outdir,norm=1)
            print(cmd)
            os.system(cmd)
    
            
