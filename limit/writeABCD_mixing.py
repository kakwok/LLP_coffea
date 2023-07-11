import numpy as np
import os
import sys
sys.path.insert(0,"../")
import pickle
import json
from optparse import OptionParser

def writeBinProcSection(text_file,signal_rate,bkg_proc):
    nBins = 4  # abcd bins
    procs = {} # name:[a,b,c,d]
    procs.update(signal_rate)
    procs.update(bkg_proc)

    binLine='bin'
    procNameLine = 'process '
    procNumLine  = 'process '
    rateNumLine  = 'rate '
    
    # one item per non-zero rate
    for i,Bin in enumerate(["chA","chB","chC","chD"]):
        for j,(procName,rates) in enumerate(procs.items()):
            i_rate = rates[i] #rate of proc j in bin i
            if i_rate>=0:
                isSignal      = -1 if (procName in signal_rate.keys()) else 1 # negative proc for signal
                binLine      +=' \t'+Bin
                procNameLine +=' \t'+procName
                procNumLine  +=' \t'+str((j+1)*(isSignal))
                if procName =="bkg":
                    rateNumLine  +=' \t 1'                                  # set rate=1 for ABCD bkg
                else:
                    rateNumLine  +=' \t {0:e}'.format(i_rate)
                
    text_file.write(binLine+"\n") 
    text_file.write(procNameLine+"\n") 
    text_file.write(procNumLine+"\n") 
    text_file.write(rateNumLine+"\n") 

def writeUncSection(text_file,signal_rate,bkg_proc,bkg_unc,sig_unc):
    ### unc = {"proc_uncX":[a,b,c,d]}

    ## match num. of proc
    procs = {} # name:[a,b,c,d]
    procs.update(signal_rate)
    procs.update(bkg_proc)

    ## bkg uncertainties
    for unc_name,unc_arr in bkg_unc.items():
        unc_text = unc_name+' \t lnN'
        if not np.any([k in unc_name for k in procs.keys()]): continue                  #skip unc that does not match any proc
        #match with all proc.
        for i,Bin in enumerate(["chA","chB","chC","chD"]):
            for j,(procName,rates) in enumerate(procs.items()):   
                i_rate = rates[i] #rate of proc j in bin i
                if i_rate>=0:                                                           #valid procs 
                    if procName in unc_name:                                            #match unc. with proc name
                        if unc_arr[i] ==0:
                            unc_text += '\t - '                                         #use - for 0 unc.
                        elif type(unc_arr[i])==type("") and len(unc_arr[i].split("/"))>1:
                            unc_text += ' \t '+unc_arr[i]                               #use string for asym unc
                        else:
                            unc_text += ' \t '+str(unc_arr[i]+1)                        #simple lnN
                    else:
                        unc_text += '\t - '                      #not this proc
        text_file.write(unc_text + ' \n')                        # write the line
    ## signal uncertainties
    for unc_name,unc_dict in sig_unc.items():
        unc_text = unc_name+' \t lnN'
        for unc_proc,unc_arr in unc_dict.items():                                           ##  unc_dict = {"HNL":[0.016,0.016, 0.016 ,0.016]}
            #print(unc_proc,unc_arr)
            for i,Bin in enumerate(["chA","chB","chC","chD"]):
                for j,(procName,rates) in enumerate(procs.items()):                          ## match proc name with unc_proc
                    i_rate = rates[i]                                                    #rate of proc j in bin i
                    if procName == unc_proc: 
                        if i_rate>=0:                                                         
                            if unc_arr[i] ==0:
                                unc_text += '\t - '                                         #use - for 0 unc.
                            elif type(unc_arr[i])==type("") and len(unc_arr[i].split("/"))>1:
                                unc_text += ' \t '+unc_arr[i]                               #use string for asym unc
                            else:
                                unc_text += ' \t '+str(unc_arr[i]+1)                        #simple lnN
                    else:
                        if i_rate>0:                                                     #skip ABC channels other bkg Zmumu_CSC
                            unc_text += '\t - '                                             #skip bkg proc
        text_file.write(unc_text + ' \n')                        # write the line
    return



def make_datacard_2sig(outDataCardsDir,modelName,  signal_rate, norm, bkg_rate, observation, bkg_unc,  sig_unc):
    ## bkg rate        = ABCD background        {bkg_name:[a,b,c,d]}
    ## signal          = ABCD signal            {sig_name:[a,b,c,d]}
    ## bkg_unc         = bkg unc               {proc_name+"_uncX":[a,b,c,d]}
    ## sig_unc         = bkg unc               {sig_name+"_unc":[a,b,c,d]}
    #a,b,c,d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    a,b,c,d = bkg_rate['bkg'][0], bkg_rate['bkg'][1], bkg_rate['bkg'][2], bkg_rate['bkg'][3]
    #c1 = a/b
    #c1 = b/a        ## DT convention
    #c2 = c/b
    #nSig = len(signal_rate.keys())
    nProc = len(signal_rate.keys())+len(bkg_rate.keys())-1
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# signal norm {0} \n'.format(norm))

    text_file.write('imax {0} \n'.format(4))
    text_file.write('jmax {0} \n'.format(nProc))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')

    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \t chD \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n'.format(observation[0],observation[1],observation[2],observation[3]))
    text_file.write('------------------------------ \n')
    writeBinProcSection(text_file,signal_rate,bkg_rate)
    text_file.write('------------------------------ \n')
    #### uncertainties ####
    writeUncSection(text_file,signal_rate,bkg_rate,bkg_unc,sig_unc)
    text_file.write('NA_val      extArg         {0}   [{0},{0}] \n'.format(a))
    text_file.write('NB_val      extArg         {0}   [{0},{0}] \n'.format(b))
    text_file.write('NC_val      extArg         {0}   [{0},{0}] \n'.format(c))
    text_file.write('NA_x        extArg          0       [-7,7] \n')
    text_file.write('NB_x        extArg          0       [-7,7] \n')
    text_file.write('NC_x        extArg          0       [-7,7] \n')
    text_file.write('NA  rateParam       chA     bkg      (@0*(1+1/sqrt(@0))**@1)       NA_val,NA_x \n')
    text_file.write('NB  rateParam       chB     bkg      (@0*(1+1/sqrt(@0))**@1)       NB_val,NB_x \n')
    text_file.write('NC  rateParam       chC     bkg      (@0*(1+1/sqrt(@0))**@1)       NC_val,NC_x \n')
    text_file.write('ND  rateParam       chD     bkg      (@0*@2/@1)                     NA,NB,NC ## D=A*C/B\n')
    text_file.write('norm  rateParam       *     HNL      1 \n')


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

def toInt(h):
    return list(h.values(True).values())[0]

def predIntimeFromOOT(h,size,dphi_lep,dphi_met,isSignal=False,kfactor=0.25,lumi=137/23.):
    # kfactor=0.25
    # lumi = 137/23.
    
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

def scale2017(h1,h2):
    ##2016+2018 = 51.2+30.9  =82.1## mu
    h1.scale(82.1)
    ##2017 = 37.9 ## mu
    h2.scale(37.9)
    h1.add(h2)      ## this is done in-place
    return 
def scale_mu(h1,h2,h3):
    ##2018 = 51.2
    #h1.scale(51.2) ##muon)
    h1.scale(59.6) ##muon)
    ##2017 = 37.0 ## muon
    h2.scale(41.5)
    ##2016 = 30.9 ## muon
    h3.scale(35.9)
    h1.add(h2)      ## this is done in-place
    h1.add(h3)      ## this is done in-place

def scale_ele(h1,h2,h3):
    ##2018 = 51.6
    h1.scale(59.6)
    ##2017 = 39.4 ## ele 
    h2.scale(41.5)
    ##2016 = 30.9 ## ele
    h3.scale(35.9)
    h1.add(h2)      ## this is done in-place
    h1.add(h3)      ## this is done in-place
    return 

def writeYields(cut = None,muon=True,outf="yields.json",debug=True,shifts=None, tauSignals=False,unblind=False,useOOT=True,mixing=True):
    from coffea import hist
    if cut==None:
        ## nHit, dphi_lep, dphi_MET
        dphi_lepcuts = np.linspace(0,np.pi,31)[1:-2]
        if muon:
            cut = {"CSC":(220,dphi_lepcuts[-10],0.7), "DT":(110,dphi_lepcuts[-8],0.7)}
        else:
            cut = {"CSC":(160,dphi_lepcuts[-2],0.7) , "DT":(100,dphi_lepcuts[-10],0.7)}
    path_prefix="/uscms_data/d3/kkwok/LLP/CMSSW_10_6_20/src/LLP_coffea/"
    if muon:
        lumi = 137
        #signals = loadSignalFromJson("./limitpoints_muon.json",True,cut,debug,lumi) ## this is too slow 
        #bkg     = loadbkg("../HNL_histograms_Feb23_muons_data.pickle",True,cut,debug)
        #signals = loadhist("../HNL_histograms_Feb23_muons_signal.pickle",True,cut,debug)
        #signals.update( loadhist("../HNL_histograms_Mar1_muons_signal.pickle",True,cut,debug))
        #bkg     = loadbkg("../HNL_histograms_Nov29_muons_all.pickle",True,cut,debug)
        #bkg      = loadbkg("../HNL_histograms_Dec6_muons_data.pickle",True,cut,debug)
        #signals  = loadhist("../HNL_histograms_Dec6_muons_signals.pickle",True,cut,debug,lumi,
        #                    "../HNL_histograms_Dec6_muons_signals_2017.pickle" )
        #bkg      = loadbkg("../HNL_histograms_Dec6_muons_noLooseID_data.pickle",True,cut,debug)
        #signals  = loadhist("../HNL_histograms_Dec6_muons_noLooseID_signals.pickle",True,cut,debug,lumi,
        #                    "../HNL_histograms_Dec6_muons_noLooseID_signals_2017.pickle" )
        #bkg      = loadbkg("../HNL_histograms_Dec14_muons_data.pickle",True,cut,debug,useOOT) ## v14
        bkg      = loadbkg(path_prefix+"HNL_histograms_Jun22_muons_data.pickle",True,cut,debug,useOOT) ## v19
        if tauSignals:
            signals  = loadhist(path_prefix+"HNL_histograms_Apr27_tau_signals_mu_2018.pickle",True,cut,debug,lumi,
                                path_prefix+"HNL_histograms_Apr27_tau_signals_mu_2017.pickle" ,
                                path_prefix+"HNL_histograms_Apr27_tau_signals_mu_2016.pickle" )
        else:
            signals  = loadhist(path_prefix+"HNL_histograms_Apr27_muons_signals_2018.pickle",True,cut,debug,lumi,
                                path_prefix+"HNL_histograms_Apr27_muons_signals_2017.pickle" ,
                                path_prefix+"HNL_histograms_Apr27_muons_signals_2016.pickle" )
    else:
        lumi = 137
        #signals = loadSignalFromJson("./limitpoints_ele.json",False,cut,debug,lumi)
        #bkg     = loadbkg("../HNL_histograms_Feb3_electrons.pickle",False,cut,debug)
        #signals = loadhist("../HNL_histograms_Apr1_ele_signal.pickle",False,cut,debug)
        #bkg     = loadbkg("../HNL_histograms_Nov29_ele_all.pickle",False,cut,debug)
        #bkg     = loadbkg("../HNL_histograms_Dec6_ele_data.pickle",False,cut,debug)
        #signals = loadhist("../HNL_histograms_Dec6_ele_signals.pickle",False,cut,debug,lumi)
        #bkg     = loadbkg("../HNL_histograms_Dec6_ele_noLooseID_data.pickle",False,cut,debug)
        #signals = loadhist("../HNL_histograms_Dec6_ele_noLooseID_signals.pickle",False,cut,debug,lumi)
        #bkg     = loadbkg("../HNL_histograms_Dec14_ele_data.pickle",False,cut,debug,useOOT)
        bkg     = loadbkg(path_prefix+"HNL_histograms_May14_ele_data.pickle",False,cut,debug,useOOT)
        if tauSignals:
            signals = loadhist(path_prefix+"HNL_histograms_Apr27_tau_signals_ele_2018.pickle",False,cut,debug,lumi,
                               path_prefix+"HNL_histograms_Apr27_tau_signals_ele_2017.pickle",
                               path_prefix+"HNL_histograms_Apr27_ele_signals_2016.pickle")
        else:
            signals = loadhist( path_prefix+"HNL_histograms_Apr27_ele_signals_2018.pickle",False,cut,debug,lumi,
                                path_prefix+"HNL_histograms_Apr27_ele_signals_2017.pickle",
                                path_prefix+"HNL_histograms_Apr27_ele_signals_2016.pickle")


    #data = {**bkg,**signals}
    ## copy data from loaded hist 
    data = {}
    for k,v in bkg.items():
        if not unblind:
            v['CSC'][-1] =-1
            v['DT'][-1]  =-1
            v['CSC_unc'][-1] =-1
            v['DT_unc'][-1]  =-1
            if "DT_MB2" in v.keys():
                v['DT_MB2'][-1]  =-1
                v['DT_MB34'][-1]  =-1
                v['DT_MB2_unc'][-1]  =-1
                v['DT_MB34_unc'][-1]  =-1
        data[k] = v 
            
    for k,v in signals.items(): data[k] = v 

    ## serialize np arrays for JSON
    for k,v in data.items():
       for region in v.keys():
            if type(v[region])==type(np.array([])):
                v.update({region:v[region].tolist()}) 

    with open(outf,"w") as f:
       f.write(json.dumps(data,indent=4))
    #if shifts is not None:
    #    shiftYields(outf,shifts)
    writeDiracYields(outf)
    
    if mixing:
        print("Mixing")
        print(muon)
        if muon:
            lumi = 137
            bkg      = loadbkg(path_prefix+"HNL_histograms_Jun22_muons_data.pickle",True,cut,debug,useOOT) ## v19 
            #signals  = loadhist(path_prefix+"HNL_histograms_Mar3_muons_signals_2018.pickle",True,cut,debug,lumi,
            #                    path_prefix+"HNL_histograms_Mar3_muons_signals_2017.pickle" ,
            #                    path_prefix+"HNL_histograms_Mar3_muons_signals_2016.pickle" )

            #signals_tau =  loadhist(path_prefix+"HNL_histograms_Mar3_tau_signals_mu_2018.pickle",True,cut,debug,lumi,
            #                    path_prefix+"HNL_histograms_Mar3_tau_signals_mu_2017.pickle" ,
            #                    path_prefix+"HNL_histograms_Mar3_tau_signals_mu_2016.pickle" )
            signals  = loadhist(path_prefix+"HNL_histograms_Apr27_muons_signals_2018.pickle",True,cut,debug,lumi,
                                path_prefix+"HNL_histograms_Apr27_muons_signals_2017.pickle" ,
                                path_prefix+"HNL_histograms_Apr27_muons_signals_2016.pickle" )
            signals_tau  = loadhist(path_prefix+"HNL_histograms_Apr27_tau_signals_mu_2018.pickle",True,cut,debug,lumi,
                                path_prefix+"HNL_histograms_Apr27_tau_signals_mu_2017.pickle" ,
                                path_prefix+"HNL_histograms_Apr27_tau_signals_mu_2016.pickle" )

        else:
            lumi=137
            bkg     = loadbkg(path_prefix+"HNL_histograms_May14_ele_data.pickle",False,cut,debug,useOOT)
            #signals = loadhist( path_prefix+"HNL_histograms_Mar3_ele_signals_2018.pickle",False,cut,debug,lumi,
            #                    path_prefix+"HNL_histograms_Mar3_ele_signals_2017.pickle",
            #                    path_prefix+"HNL_histograms_Mar3_ele_signals_2016.pickle")

            #signals_tau = loadhist(path_prefix+"HNL_histograms_Mar3_tau_signals_ele_2018.pickle",False,cut,debug,lumi,
            #                   path_prefix+"HNL_histograms_Mar3_tau_signals_ele_2017.pickle",
            #                   path_prefix+"HNL_histograms_Mar3_ele_signals_2016.pickle")

           
            signals = loadhist(path_prefix+"HNL_histograms_Apr27_tau_signals_ele_2018.pickle",False,cut,debug,lumi,
                               path_prefix+"HNL_histograms_Apr27_tau_signals_ele_2017.pickle",
                               path_prefix+"HNL_histograms_Apr27_ele_signals_2016.pickle")
        
            signals_tau = loadhist( path_prefix+"HNL_histograms_Apr27_ele_signals_2018.pickle",False,cut,debug,lumi,
                                path_prefix+"HNL_histograms_Apr27_ele_signals_2017.pickle",
                                path_prefix+"HNL_histograms_Apr27_ele_signals_2016.pickle")



        data = {}
        
        for k,v in bkg.items():
            if not unblind:
                v['CSC'][-1] =-1
                v['DT'][-1]  =-1
                v['CSC_unc'][-1] =-1
                v['DT_unc'][-1]  =-1
                if "DT_MB2" in v.keys():
                    v['DT_MB2'][-1]  =-1
                    v['DT_MB34'][-1]  =-1
                    v['DT_MB2_unc'][-1]  =-1
                    v['DT_MB34_unc'][-1]  =-1
            data[k] = v
        for signal in [signals,signals_tau]:
            for k,v in signal.items(): data[k] = v
        
        for k,v in signals.items():
            print(k)

        ## serialize np arrays for JSON
        for k,v in data.items():
            for region in v.keys():
                if type(v[region])==type(np.array([])):
                    v.update({region:v[region].tolist()})

        #print("DATA")
        #print("#"*100)
        #print(data) 
        #if shifts is not None:
        #    shiftYields("yields_mixing.json",shifts)
        writeMixingYields(data,muon,outf,shifts)
    return 

def loadbkg(fin='../HNL_histograms_Feb3_electrons.pickle',muon=False,cut=None,debug=True,useOOT=True):
    from coffea import hist

    with open(fin,'rb')  as f:
        out = pickle.load(f)
    if muon:
        datasets = [
                'Muon_2016B','Muon_2016C','Muon_2016D','Muon_2016E','Muon_2016F','Muon_2016G',"Muon_2016H",
                'Muon_2017B','Muon_2017C','Muon_2017D','Muon_2017E','Muon_2017F',
                'Muon_2018A','Muon_2018B','Muon_2018C','Muon_2018D'] 
        lumiScale = 1 
    else:
        datasets = [
            'EGamma_2016B','EGamma_2016C','EGamma_2016D','EGamma_2016E','EGamma_2016F','EGamma_2016G','EGamma_2016H',    
            'EGamma_2017B','EGamma_2017C','EGamma_2017D','EGamma_2017E','EGamma_2017F',
            'EGamma_2018A','EGamma_2018B','EGamma_2018C','EGamma_2018D']
        lumiScale = 1  

    h = out['dphi_cluster_csc'].integrate("dataset",datasets)
    hdt = out['dphi_cluster_dt'].integrate("dataset",datasets)

    if "ABCD_dt_MB2" in [s.name for s in hdt.identifiers("region")]:
            hasMB2 = True
    else:
            hasMB2 = False

    if useOOT:
        h_region    = h.integrate("region","ABCD_OOT")
        h_region_dt = hdt.integrate("region","ABCD_dt_OOT")
        if hasMB2:
            h_region_dt_MB2 = hdt.integrate("region","ABCD_dt_OOT_MB2")
            h_region_dt_MB34 = hdt.integrate("region","ABCD_dt_OOT_MB34")
        kfactorCSC=0.25
        kfactorDT=0.9
    else:
        h_region    = h.integrate("region","ABCD")
        h_region_dt = hdt.integrate("region","ABCD_dt")
        if hasMB2:
            h_region_dt_MB2 = hdt.integrate("region","ABCD_dt_MB2")
            h_region_dt_MB34 = hdt.integrate("region","ABCD_dt_MB34")
        kfactorCSC=1
        kfactorDT=1

    CSC,CSC_unc = predIntimeFromOOT(h_region,cut['CSC'][0],cut["CSC"][1],cut['CSC'][2],False,kfactorCSC,lumiScale)
    DT,DT_unc   = predIntimeFromOOT(h_region_dt,cut["DT"][0],cut["DT"][1],cut["DT"][2],False,kfactorDT,lumiScale)
    if hasMB2:
        DT_MB2,DT_MB2_unc   = predIntimeFromOOT(h_region_dt_MB2,cut["DT"][0],cut["DT"][1],cut["DT"][2],False,kfactorDT,lumiScale)
        DT_MB34,DT_MB34_unc = predIntimeFromOOT(h_region_dt_MB34,cut["DT"][0],cut["DT"][1],cut["DT"][2],False,kfactorDT,lumiScale)

    if debug:
        print("bkg CSC = " ,CSC)
        print("bkg CSCunc = " ,CSC_unc)
        print("bkg DT = " ,DT)
        print("bkg DTunc = ", DT_unc)
        if hasMB2:
                print("bkg DT_MB2 = " ,DT_MB2)
                print("bkg DTunc_MB2 = ", DT_MB2_unc)
                print("bkg DT_MB34 = " ,DT_MB34)
                print("bkg DTunc_MB34 = ", DT_MB34_unc)

    data={}
    if hasMB2:
        data["bkg"] = {"CSC":CSC,"DT":DT,"CSC_unc":CSC_unc,"DT_unc":DT_unc,
                                 "DT_MB2":DT_MB2  ,"DT_MB2_unc":DT_MB2_unc,
                                 "DT_MB34":DT_MB34,"DT_MB34_unc":DT_MB34_unc,
                                 "norm":1}
    else:
        data["bkg"] = {"CSC":CSC,"DT":DT,"CSC_unc":CSC_unc,"DT_unc":DT_unc,"norm":1}
    return data 

def scaleHist(out,signal_name,lumi,cut,debug):
    from coffea import hist
    
    for k,h in out.items():
        if (type(h)!=hist.Hist): continue
        h.scale({ d: lumi for d in h.identifiers("dataset")} , axis="dataset") 

    h = out['dphi_cluster_csc']
    signal = h.integrate("dataset",signal_name).integrate("region","ABCD")
    hdt = out['dphi_cluster_dt']
    signal_dt = hdt.integrate("dataset",signal_name).integrate("region","ABCD_dt")

    ## not used for signal
    kfactor=1
    lumi=1
    CSC,CSC_unc = predIntimeFromOOT(signal,cut["CSC"][0],cut["CSC"][1],cut["CSC"][2],True,kfactor,lumi)
    DT,DT_unc   = predIntimeFromOOT(signal_dt,cut["DT"][0],cut["DT"][1],cut["DT"][2],True,kfactor,lumi)
    if "rwctau" in signal_name:
        ct = signal_name.split("_")[-1].replace("rwctau","pl")
        sample_name = ("_".join(signal_name.split("_")[:-2]+[ct]))
    else:
        sample_name = signal_name 
    if debug:
        print(signal_name," CSC = " ,CSC)
        print(signal_name," CSCunc = " ,CSC_unc)
        print(signal_name," DT = " ,DT)
        print(signal_name," DTunc = ", DT_unc)
    return {"CSC":CSC,"DT":DT,"CSC_unc":CSC_unc,"DT_unc":DT_unc,"norm":1}
 

### load signals from json
def loadSignalFromJson(fin='../signals_Nov18_ele.json',muon=True,cut=None,debug=True,lumi=137):
    from coffea.nanoevents import NanoEventsFactory, NanoAODSchema,BaseSchema
    from HNLprocessor.HNLproc_4 import MyProcessor

    import time
    proc    = MyProcessor(not(muon))
    with open(fin,'r') as f:
        signals = json.load(f)

    data = {} 

    for dataset,fpath in signals.items():
        start = time.perf_counter()
        if len(fpath)>1:
            print("Skipping %s, because it has more than 1 file"%dataset)
        else:
            if debug:
                print("dataset = ", dataset)
                print("     Loading events  ...")
                tic = time.perf_counter()
            events = NanoEventsFactory.from_root(fpath[0], schemaclass=BaseSchema,treepath='MuonSystem', metadata={"dataset":dataset }).events()
            if debug: print("     Processing events  ...")
            out = proc.process(events)
            out = proc.postprocess(out)
            if debug:
                toc = time.perf_counter()
                #print(f"     processing time =   {toc - tic:0.4f} seconds")
            del(events)
            data[dataset] = scaleHist(out,dataset,lumi,cut,debug) 
    stop = time.perf_counter()
    #if debug: print(f"     total processing time =   {stop - start:0.4f} seconds")
    return data

    
def loadhist(fin='../HNL_histograms_Feb23_muons_signal.pickle',muon=True,cut=None,debug=True,lumi=137,fin2017=None,fin2016=None):
    from coffea import hist
    
    with open(fin,'rb')  as f:
        out = pickle.load(f)
    
    if fin2017 is not None and fin2016 is not None:
        with open(fin2017,'rb')  as f:
            out_2016 = pickle.load(f)
        with open(fin2017,'rb')  as f:
            out_2017 = pickle.load(f)
        for k,h in out.items():
            if (type(h)!=hist.Hist): continue
            h_2017 = out_2017[k]    # find the same histogram in 2017 output
            h_2016 = out_2016[k]    # find the same histogram in 2016 output
            if muon:
                scale_mu(h,h_2017,h_2016)
            else:
                scale_ele(h,h_2017,h_2016)
    elif fin2017 is not None:
        with open(fin2017,'rb')  as f:
            out_2017 = pickle.load(f)
        for k,h in out.items():
            if (type(h)!=hist.Hist): continue
            h_2017 = out_2017[k]    # find the same histogram in 2017 output
            scale2017(h,h_2017)
    else:
        for k,h in out.items():
            if (type(h)!=hist.Hist): continue
            #print("scaling ",k)
            h.scale({ d: lumi for d in h.identifiers("dataset")} , axis="dataset") 

    signalNames = [ s.name for s in out['dphi_cluster_csc'].identifiers("dataset")]   
    data = {} 
    for signal_name in signalNames:
        h = out['dphi_cluster_csc']
        signal = h.integrate("dataset",signal_name).integrate("region","ABCD")
        hdt = out['dphi_cluster_dt']
        signal_dt = hdt.integrate("dataset",signal_name).integrate("region","ABCD_dt")
        if "ABCD_dt_MB2" in [s.name for s in hdt.identifiers("region")]:
            hasMB2 = True
        else:
            hasMB2 = False
        if hasMB2:
            signal_dt_MB2 = hdt.integrate("dataset",signal_name).integrate("region","ABCD_dt_MB2")
            signal_dt_MB34 = hdt.integrate("dataset",signal_name).integrate("region","ABCD_dt_MB34")

        ### NOT USED FOR SIGNAL
        lumi = 1
        kfactor=1     ## Muon, CSC InT/OOT 

        ## CSC cuts
        #CSC_cls_SF = (1-0.015)*(1-0.033)*(1-0.060)   ## ~0.9
        CSC_cls_SF = (1-0.021)*(1-0.028)*(1-0.068)   ## 0.88 
        CSC,CSC_unc = predIntimeFromOOT(signal,cut["CSC"][0],cut["CSC"][1],cut["CSC"][2],True,kfactor,lumi)
        CSC     *=  CSC_cls_SF
        CSC_unc *=  CSC_cls_SF

        ### DT cuts
        DT_cls_SF = 0.907*0.932*0.910        ## ~0.76
        DT_cls_SF = 0.932*0.910
        DT,DT_unc = predIntimeFromOOT(signal_dt,cut["DT"][0],cut["DT"][1],cut["DT"][2],True,kfactor,lumi)
        DT *=DT_cls_SF
        DT_unc*=DT_cls_SF
        if hasMB2:
            DT_MB2,DT_MB2_unc = predIntimeFromOOT(signal_dt_MB2,cut["DT"][0],cut["DT"][1],cut["DT"][2],True,kfactor,lumi)
            DT_MB34,DT_MB34_unc = predIntimeFromOOT(signal_dt_MB34,cut["DT"][0],cut["DT"][1],cut["DT"][2],True,kfactor,lumi)
            DT_MB2 *= DT_cls_SF
            DT_MB2_unc *= DT_cls_SF
            DT_MB34 *= DT_cls_SF
            DT_MB34_unc *= DT_cls_SF
        if "rwctau" in signal_name:
            ct = signal_name.split("_")[-1].replace("rwctau","pl")
            sample_name = ("_".join(signal_name.split("_")[:-2]+[ct]))
        else:
            sample_name = signal_name 
        #if "1p0" in sample_name:
        if debug:
            print(signal_name," CSC = " ,CSC)
            print(signal_name," CSCunc = " ,CSC_unc)
            print(signal_name," DT = " ,DT)
            print(signal_name," DTunc = ", DT_unc)
            if hasMB2:
                print(signal_name," DT_MB2 = " ,DT_MB2)
                print(signal_name," DTunc_MB2 = ", DT_MB2_unc)
                print(signal_name," DT_MB34 = " ,DT_MB34)
                print(signal_name," DTunc_MB34 = ", DT_MB34_unc)
        #data[sample_name] = {"CSC":CSC,"DT":DT,"norm":1}
        if hasMB2:
            data[sample_name] = {"CSC":CSC,"DT":DT,"CSC_unc":CSC_unc,"DT_unc":DT_unc,
                                 "DT_MB2":DT_MB2  ,"DT_MB2_unc":DT_MB2_unc,
                                 "DT_MB34":DT_MB34,"DT_MB34_unc":DT_MB34_unc,
                                 "norm":1}
        else:
            data[sample_name] = {"CSC":CSC,"DT":DT,"CSC_unc":CSC_unc,"DT_unc":DT_unc,"norm":1}
    return data 

def makeAllcards(f_yield,outdir="./combine/HNL_datacards/",suffix="",test=False):

    from collections import OrderedDict
    if not os.path.exists(outdir):
        print("mk dir = ", outdir)
        os.mkdir(outdir)
    norm = 100.0    # 1% BR
    suffix = ""
    
    with open(f_yield,'r') as f:
        data = json.load(f)
        if "DT_MB2" in data['bkg'].keys():
            hasMB2=True
        else:
            hasMB2 = False
    

    #ZmumuCR = CSC:73      DT:172  (200,130)
    #ZmumuCR = CSC:54      DT:130  (220,150)

    if options.muon:
        bkg_proc_CSC= {
         #   "Zmumu_CSC":[-1,-1,-1,2.07], ## 2.5% TF 
            #"Zmumu_CSC":[-1,-1,-1,4.68], ## 5.6% TF 
            #"Zmumu_CSC":[-1,-1,-1,4.15], ## 5.69% TF, 200 cut
        }
        #bkg_proc_DT= {
        # #   "Zmumu_DT":[-1,-1,-1,1.95],    # 1.0% TF
        #    #"Zmumu_DT":[-1,-1,-1,2.56],     # 1.3% TF
        #    #"Zmumu_DT":[-1,-1,-1,1.69],     # 1.3% TF, 150
        #    #"Zmumu_DT":[-1,-1,-1,5.69],     # 1.7% TF, 130, loose ID
        #    "Zmumu_DT":[-1,-1,-1,2.45],     # 1.47% TF, 130, no loose ID
        #    #"Zmumu_DT":[-1,-1,-1,23.94],     # 14.3% high nHit TF, 130, no loose ID
        #}
        ## must use convention :     procName_uncName
        #bkg_unc = {
        #    #"bkg_sys_D":[0  ,0 ,0 ,0.5],
        #    #"Zmumu_CSC_sys_D":[0 ,0 ,0 ,0.32],      
        #    #"Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.2 ],
        #    "Zmumu_CSC_sys_D":[0 ,0 ,0 ,0.29], ## 5.69% TF, 200 cut 
        #    #"Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.24], # 1.7% TF, 130, loose ID
        #    "Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.3], # 1.47% TF, 130, no loose ID
        #    #"Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.75], # 14.3% high nHit TF 
        #    #"Zmumu_DT_MB2_sys_D" :[0 ,0 ,0 ,0.42],     # 0.42=3/7
        #}
       #bkg_proc_CSC= {    "Zmumu_CSC":[-1,-1,-1,3.5]} ## 4.7% TF, 200 cut}
        bkg_proc_CSC= {    "Zmumu_CSC":[-1,-1,-1,4.1]} ## 4.8% TF, 200 cut, full lumi}
        #bkg_unc_CSC = {    "Zmumu_CSC":[-1,-1,-1,0.26]    }
        bkg_unc_CSC = {    "Zmumu_CSC":[-1,-1,-1,0.30]    } #new systematics
        ## single-bin DT with total Zmumu
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
        "dt_mixingSyst"   :{       "HNL":[0.10,0.10, 0.10 ,0.10],    },
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

    with open(f_yield,'r') as f:
        data = json.load(f)
        bkg_rate_CSC = {}       # {bkg: [a,b,c,d]}
        bkg_rate_CSC.update({"bkg":data['bkg']["CSC"]})          # {bkg: [a,b,c,d]}
        bkg_rate_CSC.update(bkg_proc_CSC)                        # add other bkg 
        bkg_rate_DT = {}       # {bkg: [a,b,c,d]}
        bkg_rate_DT.update({"bkg":data['bkg']["DT"]})           # {bkg: [a,b,c,d]}
        bkg_rate_DT.update(bkg_proc_DT)                         # add other bkg 
        if hasMB2:
            bkg_rate_DT_MB2 = {} 
            bkg_rate_DT_MB2.update({"bkg":data['bkg']["DT_MB2"]})
            bkg_rate_DT_MB2.update(bkg_proc_DT_MB2)
            bkg_rate_DT_MB34 = {} 
            bkg_rate_DT_MB34.update({"bkg":data['bkg']["DT_MB34"]})
            bkg_rate_DT_MB34.update(bkg_proc_DT_MB34)
            
        for name,signal in data.items():
            if "bkg" in name: continue
            #if not ("2p0" in name or "2p1" in name):  continue
            if not ("pl1000" in name): continue
            if suffix:
                name = name+suffix
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
            if hasMB2:
                obs = bkg_rate_DT_MB2['bkg']
                sigRate = {"HNL":np.array(signal["DT_MB2"]) /norm}
                #obs[-1] = bkg_rate_DT_MB2['bkg'][0]*bkg_rate_DT_MB2['bkg'][2]/bkg_rate_DT_MB2['bkg'][1]  ## force D=A*C/B
                make_datacard_2sig(outdir,name+"_DT_MB2", sigRate, norm, bkg_rate_DT_MB2, obs, bkg_unc_DT_MB2,  sig_unc_DT)
                obs = bkg_rate_DT_MB34['bkg']
                sigRate = {"HNL":np.array(signal["DT_MB34"]) /norm}
                #obs[-1] = bkg_rate_DT_MB34['bkg'][0]*bkg_rate_DT_MB34['bkg'][2]/bkg_rate_DT_MB34['bkg'][1]  ## force D=A*C/B
                make_datacard_2sig(outdir,name+"_DT_MB34", sigRate, norm, bkg_rate_DT_MB34, obs, bkg_unc_DT_MB34,  sig_unc_DT)
   
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
                    
            RunLimit(outdir,name,"CSC",None,None,test)
            if hasMB2:
                #RunLimit(outdir,name,"DT_MB2",None,None,test)
                #RunLimit(outdir,name,"DT_MB34",None,None,test)
                RunLimit(outdir,name,"DT_MB2","DT_MB34","DTcomb",test)
                pass
            else:
                RunLimit(outdir,name,"DT",None,None,test)

            #Combination
            if hasMB2:
                RunLimit(outdir,name,"CSC","DTcomb","comb",test)       # CSC+MB2+MB34
                #RunLimit(outdir,name,"CSC","DT_MB34","combNoMB2",test) # CSC+MB34
                #RunLimit(outdir,name,"CSC","DT","combSingleDTBin",test)            # CSC+DT
                #RunLimit(outdir,name,"CSC","DT_MB34","comb",test) # CSC+MB34
                pass
            else:
                RunLimit(outdir,name,"CSC","DT","comb",test)
            
    return 

               
import sys
sys.path.insert(0,"../")

# Find the corresponding yield of (m,ct)-> (m,xx)
def shift_ctau(N_yield,m_old,ct_old,m_new,forTauHNL=False):
    from HNLprocessor.util import f_1m,f_xsec, f_1m_tau,f_xsec_tau 
    ct_new = (m_new/m_old) * ct_old
    if forTauHNL:
        return N_yield * f_xsec_tau(m_new)(ct_new)/f_xsec_tau(m_old)(ct_old)
    else:
        return N_yield * f_xsec(m_new)(ct_new)/f_xsec(m_old)(ct_old)

def shift_ctau_mixing(N_yield,m_old,ct_old,m_new,name):
    from HNLprocessor.util import f_1m,f_xsec, f_1m_tau,f_xsec_tau
    ct_new = (m_new/m_old) * ct_old
    fe = float(name.split("-")[1][2:].replace("p","."))
    fmu = float(name.split("-")[2][3:].replace("p","."))
    ftau = float(name.split("-")[3].split("_")[0][4:].replace("p","."))
    return (ftau*(N_yield * f_xsec_tau(m_new)(ct_new)/f_xsec_tau(m_old)(ct_old)))+((fe+fmu)*(N_yield * f_xsec(m_new)(ct_new)/f_xsec(m_old)(ct_old)))  
  

# Dirac yield = 2x Majorana yield 
def writeDiracYields(f_yield,shifts=None):
    dirac={}
    with open(f_yield,'r') as f:
        data = json.load(f)
        for name,signal in data.items():
            if name=="bkg":
                dirac['bkg'] = signal # copy the bkg
            else:
                diracName = name.replace("Type","DiracType")
                # double the yields and unc.
                dirac[diracName] = {region:(np.array(yields)*2).tolist() for region,yields in signal.items()}
    with open(f_yield.replace(".json","_dirac.json"),'w') as f:
        f.write(json.dumps(dirac,indent=4))             
    if shifts is not None:
            print("Shifiting yields for mixed")
            shiftYields(f_yield.replace(".json","_dirac.json"),shifts)
    return 

def writeMixingYields(data,muon,f_yield,shifts):
    mixing = {}
    #mixing_points = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    mixing_points = np.arange(0,1.02,0.02)
    for fe in mixing_points:
        
        for fmu in mixing_points:
      
            
      
            for ftau in mixing_points:
               
                sum_of_floats = fe+fmu+ftau
                sum_of_floats = round(sum_of_floats,2)
                
                if sum_of_floats!=1: continue
                
               
                for name,signal in data.items():
                        #if "mHNL2p0" not in name: continue
                        if "bkg" in name:
                            mixing[name] = signal # copy the bkg
                        else:
                            mixed_sample_name = "HNL_mixed-fe{}-fmu{}-ftau{}_".format(round(fe,2),round(fmu,2),round(ftau,2)).replace(".","p")+name[name.find("mHNL"):]
                            print(mixed_sample_name)                    
                            if muon:
                                if "muon" in name:
                                    muon_mass_lifetime = name[name.find("mHNL")+4:]
                                    flag_found_match = False
                                    for tau_name,tau_signal in data.items(): 
                                        
                                        mass_lifetime = tau_name[tau_name.find("mHNL")+4:] 
                                        if "tau" in tau_name and muon_mass_lifetime==mass_lifetime:
                                            flag_found_match = True
                                            #if np.sum(np.array(tau_signal["CSC"])*ftau+np.array(signal["CSC"])*fmu)==0: continue
                                            #if np.sum(np.array(tau_signal["DT"])*ftau+np.array(signal["DT"])*fmu)==0: continue

                                            mixing[mixed_sample_name] = {muon[0]:(np.array(muon[1])*fmu+np.array(tau[1])*ftau).tolist() for muon,tau in zip(signal.items(),tau_signal.items())}
                                    
                                    if not flag_found_match:
                                        #if np.sum(np.array(signal["CSC"])*fmu)==0: continue
                                        #if np.sum(np.array(signal["DT"])*fmu)==0: continue
                                        mixing[mixed_sample_name] = {region:(np.array(yields)*fmu).tolist() for region,yields in signal.items()}
                            else:
                                #print(name)                       
                                if "electron" in name:
                                    
                                    ele_mass_lifetime = name[name.find("mHNL")+4:]
                                    print(ele_mass_lifetime)
                                    #print(signal.items())
                                    flag_found_match = False
                                    for tau_name,tau_signal in data.items():
                                        
                                        mass_lifetime = tau_name[tau_name.find("mHNL")+4:]
                                       
                                        if "tau" in tau_name and ele_mass_lifetime==mass_lifetime:
                                            flag_found_match = True
                                            print("Match found")
                                            if np.sum(np.array(tau_signal["CSC"])*ftau+np.array(signal["CSC"])*fe)==0:
                                                print("CSC yield is 0")
                                                #continue 
                                            if np.sum(np.array(tau_signal["DT"])*ftau+np.array(signal["DT"])*fe)==0:
                                                print("DT yield is 0")
                                                #continue 
                                            
                                            
                                            mixing[mixed_sample_name] = {ele[0]:(np.array(ele[1])*fe+np.array(tau[1])*ftau).tolist() for ele,tau in zip(signal.items(),tau_signal.items())}
                                    if not flag_found_match:
                                        print("No match found")
                                        if np.sum(np.array(signal["CSC"])*fe)==0: 
                                            print("CSC yield is 0")
                                            #continue  
                                        if np.sum(np.array(signal["DT"])*fe)==0: 
                                            print("DT yield is 0")
                                            #continue 
                                        mixing[mixed_sample_name] = {region:(np.array(yields)*fe).tolist() for region,yields in signal.items()}
    #if muon:
    #    with open(f_yield.replace(".json","_mixed.json"),'w') as f: 
    #        f.write(json.dumps(mixing,indent=4))
    #        print("Wrote out: ",f_yield.replace(".json","_mixed.json"))
    #else:
    #    with open(f_yield.replace(".json","_mixed.json"),'w') as f:
    #        f.write(json.dumps(mixing,indent=4))
    #        print("Wrote out: ",f_yield.replace(".json","_mixed.json"))
    with open(f_yield.replace(".json","_mixed.json"),'w') as f:
            f.write(json.dumps(mixing,indent=4))
            print("Wrote out: ",f_yield.replace(".json","_mixed.json"))
    if shifts is not None:
            print("Shifiting yields for mixed")
            shiftYields(f_yield.replace(".json","_mixed.json"),shifts)


    writeDiracYields(f_yield.replace(".json","_mixed.json"),shifts)
    
    return

# shift yielsd in f_yield.json according to shifts[{"m_src":int,"m_target":int}]
def shiftYields(f_yield,shifts):
    ## init
    print("Shifting")
    print(shifts)
    data_shifts = {}
    for shift in shifts:
        shift['CSC_target']=[]
        shift['CSC_unc_target']=[]
        shift['DT_target']=[]
        shift['DT_unc_target']=[]
        shift['DT_MB2_target']=[]
        shift['DT_MB2_unc_target']=[]
        shift['DT_MB34_target']=[]
        shift['DT_MB34_unc_target']=[]
        shift['names']=[]
    with open(f_yield,'r') as f:
        data = json.load(f)
        for name,signal in data.items():
            if name=="bkg": 
                data_shifts[name]=signal
                continue
            if "DT_MB2" in signal.keys():
                hasMB2=True
            else:
                hasMB2=False
            m =float(name.split("_")[-2].replace("mHNL","").replace("p","."))
            ct =float(name.split("_")[-1].replace("pl",""))    
            for shift in shifts:
                if m==shift["m_src"]:     
                    ct_new = (shift['m_target']/m) * ct            
                    shift["CSC_target"].append(shift_ctau(np.array(signal['CSC']),m,ct,shift['m_target'],shift['isTau']))   
                    shift["CSC_unc_target"].append(shift_ctau(np.array(signal['CSC_unc']),m,ct,shift['m_target'],shift['isTau']))            
                    shift["DT_target"].append(shift_ctau(np.array(signal['DT']),m,ct,shift['m_target'],shift['isTau']))                        
                    shift["DT_unc_target"].append(shift_ctau(np.array(signal['DT_unc']),m,ct,shift['m_target'],shift['isTau']))                      
                    if hasMB2: 
                        shift["DT_MB2_target"].append(    shift_ctau(np.array(signal['DT_MB2']),m,ct,shift['m_target'],shift['isTau']))                        
                        shift["DT_MB2_unc_target"].append(shift_ctau(np.array(signal['DT_MB2_unc']),m,ct,shift['m_target'],shift['isTau']))                        
                        shift["DT_MB34_target"].append(    shift_ctau(np.array(signal['DT_MB34']),m,ct,shift['m_target'],shift['isTau']))                        
                        shift["DT_MB34_unc_target"].append(shift_ctau(np.array(signal['DT_MB34_unc']),m,ct,shift['m_target'],shift['isTau']))                        
                    m_src = str(shift["m_src"]).replace(".","p")
                    m_tar = str(shift["m_target"]).replace(".","p") 
                    ct_src = name.split("_")[-1]
                    ct_tar = "pl"+str(int(ct_new))
                    shift["names"].append(name.replace(m_src,m_tar).replace(ct_src,ct_tar))
    for shift in shifts:
        if len(shift['names'])>0:
            for i,name in enumerate(shift['names']):
                print(name)
                data[name]={
                    "CSC":shift["CSC_target"][i].tolist(),
                    "CSC_unc":shift["CSC_unc_target"][i].tolist(),
                    "DT":shift["DT_target"][i].tolist(),  
                    "DT_unc":shift["DT_unc_target"][i].tolist(),  
                    "norm":1,
                } 
                data_shifts[name]={
                    "CSC":shift["CSC_target"][i].tolist(),
                    "CSC_unc":shift["CSC_unc_target"][i].tolist(),
                    "DT":shift["DT_target"][i].tolist(),
                    "DT_unc":shift["DT_unc_target"][i].tolist(),
                    "norm":1,
                }
                if hasMB2:
                    data[name].update( {"DT_MB2":shift["DT_MB2_target"][i].tolist()})
                    data[name].update( {"DT_MB34":shift["DT_MB34_target"][i].tolist()})
                    data_shifts[name].update( {"DT_MB2":shift["DT_MB2_target"][i].tolist()})
                    data_shifts[name].update( {"DT_MB34":shift["DT_MB34_target"][i].tolist()})
        else:
            print("Cannot find source data of mass = ",shift['m_src'])
    print(f_yield)
    with open(f_yield.replace(".json","_shifted.json"),'w') as f:
        print("write shifts out")
        #print(data)
        f.write(json.dumps(data_shifts,indent=4)) 
    return 


def shiftYields_mixed(f_yield,shifts):
    ## init
    for shift in shifts:
        shift['CSC_target']=[]
        shift['CSC_unc_target']=[]
        shift['DT_target']=[]
        shift['DT_unc_target']=[]
        shift['DT_MB2_target']=[]
        shift['DT_MB2_unc_target']=[]
        shift['DT_MB34_target']=[]
        shift['DT_MB34_unc_target']=[]
        shift['names']=[]
    with open(f_yield,'r') as f:
        data = json.load(f)
        for name,signal in data.items():
            if name=="bkg": continue
            if "DT_MB2" in signal.keys():
                hasMB2=True
            else:
                hasMB2=False
            m =float(name.split("_")[-2].replace("mHNL","").replace("p","."))
            ct =float(name.split("_")[-1].replace("pl",""))    
            for shift in shifts:
                if m==shift["m_src"]:     
                    ct_new = (shift['m_target']/m) * ct            
                    shift["CSC_target"].append(shift_ctau_mixing(np.array(signal['CSC']),m,ct,shift['m_target'],name))   
                    shift["CSC_unc_target"].append(shift_ctau_mixing(np.array(signal['CSC_unc']),m,ct,shift['m_target'],name))            
                    shift["DT_target"].append(shift_ctau_mixing(np.array(signal['DT']),m,ct,shift['m_target'],name))                        
                    shift["DT_unc_target"].append(shift_ctau_mixing(np.array(signal['DT_unc']),m,ct,shift['m_target'],name))                      
                    if hasMB2: 
                        shift["DT_MB2_target"].append(    shift_ctau_mixing(np.array(signal['DT_MB2']),m,ct,shift['m_target'],name))                        
                        shift["DT_MB2_unc_target"].append(shift_ctau_mixing(np.array(signal['DT_MB2_unc']),m,ct,shift['m_target'],name))                        
                        shift["DT_MB34_target"].append(    shift_ctau_mixing(np.array(signal['DT_MB34']),m,ct,shift['m_target'],name))                        
                        shift["DT_MB34_unc_target"].append(shift_ctau_mixing(np.array(signal['DT_MB34_unc']),m,ct,shift['m_target'],name))                        
                    m_src = str(shift["m_src"]).replace(".","p")
                    m_tar = str(shift["m_target"]).replace(".","p") 
                    ct_src = name.split("_")[-1]
                    ct_tar = "pl"+str(int(ct_new))
                    shift["names"].append(name.replace(m_src,m_tar).replace(ct_src,ct_tar))
    for shift in shifts:
        if len(shift['names'])>0:
            for i,name in enumerate(shift['names']):
                data[name]={
                    "CSC":shift["CSC_target"][i].tolist(),
                    "CSC_unc":shift["CSC_unc_target"][i].tolist(),
                    "DT":shift["DT_target"][i].tolist(),  
                    "DT_unc":shift["DT_unc_target"][i].tolist(),  
                    "norm":1,
                }  
                if hasMB2:
                    data[name].update( {"DT_MB2":shift["DT_MB2_target"][i].tolist()})
                    data[name].update( {"DT_MB34":shift["DT_MB34_target"][i].tolist()})
        else:
            print("Cannot find source data of mass = ",shift['m_src'])

    with open(f_yield,'w') as f:
        f.write(json.dumps(data,indent=4)) 
    return 



if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--loadhist', dest='loadhist', action='store_true',default = False, help='Load result from pickle')
    parser.add_option('-f', dest='fin', default = "", help='pickle name ')
    parser.add_option('--test', dest='test', action='store_true',default = False, help='test writing cards')
    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='dryRun')
    parser.add_option('--loadbkg', dest='loadbkg', action='store_true',default = False, help='Load bkg from pickle')
    parser.add_option('--writeYields', dest='writeYields', action='store_true',default = False, help='write yield fields')
    parser.add_option('--useOOT', dest='useOOT', action='store_true',default = False, help='useOOT for bkg ABCD prediction, use real ABC otherwise ')
    parser.add_option('--unblind', dest='unblind', action='store_true',default = False, help='do not blind bin D if true ')
    parser.add_option('--muon', dest='muon', action='store_true',default = False, help='make muon datacard')
    parser.add_option('--tauSignals', dest='tauSignals', action='store_true',default = False, help='make tau signals datacard')
    (options, args) = parser.parse_args()

    dphi_lepcuts = np.linspace(0,np.pi,31)[1:-2]
    #[0.10471976 0.20943951 0.31415927 0.41887902 0.52359878 0.62831853
    # 0.73303829 0.83775804 0.9424778  1.04719755 1.15191731 1.25663706
    # 1.36135682 1.46607657 1.57079633 1.67551608 1.78023584 1.88495559
    # 1.98967535 2.0943951  2.19911486 2.30383461 2.40855437 2.51327412
    # 2.61799388 2.72271363 2.82743339 2.93215314]
    ## default cuts
    if options.muon:
        cut = {"CSC":(220,dphi_lepcuts[-10],0.7), "DT":(110,dphi_lepcuts[-8],0.7)}
    else:
        cut = {"CSC":(160,dphi_lepcuts[-2],0.7), "DT":(100,dphi_lepcuts[-10],0.7)}

    #writeYields(cut,options.muon) 
    if options.loadhist:
        data = loadhist(options.fin,options.muon)
        for k,v in data.items():
            print(k,v)
    elif options.loadbkg:
        data = loadbkg()
    elif options.test:
        print("Testing script ")
        outdir = "./test/"   ### 
        cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        isMuon=options.muon
        f_yield = "yields.json"
        #f_yield = "Mixed_yields_muon.json"
        #makeAllcards(f_yield,outdir,"",True)
        shifts={}
        if options.writeYields: writeYields(cut,isMuon,f_yield,True,shifts,options.tauSignals,options.unblind,options.useOOT) 
    else:
        ### electron chan using muon cuts
        #outdir = "./combine/HNL_datacards/ele_muCuts/"
        #cut = {"CSC":(220,dphi_lepcuts[-10],0.7), "DT":(110,dphi_lepcuts[-8],0.7)}
        #outdir = "./combine/HNL_datacards/ele_v1//"
        #cut = {"CSC":(160,dphi_lepcuts[-2],0.7), "DT":(100,dphi_lepcuts[-10],0.7)}
        #########################################
        #outdir = "./combine/HNL_datacards/ele_v2/"
        #cut = {"CSC":(140,dphi_lepcuts[5],0.9), "DT":(80,dphi_lepcuts[5],0.6)}
        #isMuon=False
        #########################################
        #outdir = "./combine/HNL_datacards/ele_v3_fixed/"
        #cut = {"CSC":(140,1.0,0.7), "DT":(100,1.0,0.7)}
        #isMuon=False
        #########################################
        #print("Working on electron Channel: ")
        #outdir = "./combine/HNL_datacards/ele_v5/"   ###  v5 == v4 + 50% unc + toysFreq
        #cut = {"CSC":(160,dphi_lepcuts[-4],None), "DT":(100,dphi_lepcuts[-4],None)}
        #isMuon=False
        #########################################
        #outdir = "./combine/HNL_datacards/ele_v6/"   ###  v6 == v5 + 50% unc + toysFreq
        #cut = {"CSC":(200,dphi_lepcuts[-4],None), "DT":(130,dphi_lepcuts[-4],None)}
        #########################################
        #outdir = "./combine/HNL_datacards/ele_v7/"   ###  v7 == v6 with dphilep=2.82
        #cut = {"CSC":(200,dphi_lepcuts[-2],None), "DT":(130,dphi_lepcuts[-2],None)}
        #########################################
        print("Working on electron Channel: ")
        #outdir = "./combine/HNL_datacards/ele_v8/"   ###  v8, full run 2, looseID, new timing,
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/ele_v9/"   ###  v9, full run 2,no looseID, new timing
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/ele_v10/"   ###  v10, new xSec, new ele pT cuts, trig. bits, trig SF
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/tau_v1/ele/"   ### tau-v1 
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/ele_v11/"   ###  v11, new 2018 ele pT cuts, ID SF fixes (Jan4 pickles)
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/ele_v12/"   ###  v12, 150DT cut (Jan4 pickles) OOTxTF
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #outdir = "./combine/HNL_datacards/ele_v13/"   ###  v13, 150DT cut (Jan4 pickles) ABCD, CSC,DT cls SF 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        outdir = "./combine/HNL_datacards/ele_v20_full_lumi/"   ###  v14, 150DT cut (Jan4 pickles) ABCD, unblind 
        cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #########################################
        #outdir = "./combine/HNL_datacards/tau_v2/ele/"   ### tau v2: 150 DT, real ABC
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #outdir = "./combine/HNL_datacards/tau_v3/ele/"   ### tau v3 (as in v2): CSC,DT cls SF 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #outdir = "./combine/HNL_datacards/tau_v4/ele/"   ### tau v4 (as in v3): unblinded
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        isMuon=False

        f_yield = outdir+"yields.json"
       
        print("        yields = ", f_yield)
        if options.tauSignals:
            # shifts for tau HNLs
            shifts = [
                #{"m_src":2.0,"m_target":2.5,"isTau":True},
                {"m_src":2.0,"m_target":2.2,"isTau":True},
                {"m_src":2.0,"m_target":2.1,"isTau":True},
                {"m_src":2.0,"m_target":1.8,"isTau":True},
                {"m_src":2.0,"m_target":1.7,"isTau":True},
                {"m_src":2.0,"m_target":1.6,"isTau":True},
                {"m_src":2.0,"m_target":1.5,"isTau":True},
                {"m_src":1.0,"m_target":1.3,"isTau":True},
            ]
        else:
            shifts = [
                {"m_src":4.0,"m_target":2.5,"isTau":False},
                {"m_src":4.0,"m_target":2.8,"isTau":False},
                {"m_src":4.0,"m_target":3.0,"isTau":False},
                {"m_src":4.0,"m_target":3.1,"isTau":False},
                {"m_src":4.0,"m_target":3.2,"isTau":False},
                {"m_src":4.0,"m_target":3.3,"isTau":False},
                {"m_src":4.0,"m_target":3.4,"isTau":False},
                {"m_src":4.0,"m_target":3.5,"isTau":False},
            ]
        shifts = [
                #{"m_src":4.0,"m_target":2.5,"isTau":True},
                {"m_src":2.0,"m_target":1.5,"isTau":True},
                ]
        if not options.muon:
            print(shifts)
            if options.writeYields: writeYields(cut,isMuon,f_yield,True,shifts,options.tauSignals,options.unblind,options.useOOT) 
            else:
                makeAllcards(f_yield.replace(".json","_mixed.json"),outdir,"",options.dryRun)
                #makeAllcards(f_yield.replace(".json","_dirac.json"),outdir,"",options.dryRun)
        #########################################
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v3/"
        #cut = {"CSC":(220,dphi_lepcuts[9],1.0), "DT":(110,dphi_lepcuts[9],1.0)}
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v2/"
        #cut = {"CSC":(250,dphi_lepcuts[9],1.0), "DT":(140,dphi_lepcuts[9],1.0)}
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v1/"
        #cut = {"CSC":(220,dphi_lepcuts[-10],0.7), "DT":(110,dphi_lepcuts[-8],0.7)}
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v5/"   ### v5 == v4 + 50% unc + toysFreq
        #cut = {"CSC":(200,dphi_lepcuts[-4],None), "DT":(120,dphi_lepcuts[-4],None)}
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v6/"   ### v6 == v4 + 50% unc in D+ toysFreq
        #cut = {"CSC":(200,dphi_lepcuts[-4],None), "DT":(130,dphi_lepcuts[-4],None)}
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v7/"   ### v7 == v5 + but with 2.82
        #cut = {"CSC":(200,dphi_lepcuts[-2],None), "DT":(130,dphi_lepcuts[-2],None)}
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v8/"   ### v8 == v7 + Zmumubkg (TF 2.5%/1%) 
        #outdir = "./combine/HNL_datacards/muon_v9/"   ### v9 == v7 + Zmumubkg (TF 5.6%/1.3%) 
        #cut = {"CSC":(200,dphi_lepcuts[-2],None), "DT":(130,dphi_lepcuts[-2],None)}
        #outdir = "./combine/HNL_datacards/muon_v10/"   ### v10 == v9 + Zmumubkg (TF 5.6%/1.3%) at 150/220 
        #cut = {"CSC":(220,dphi_lepcuts[-2],None), "DT":(150,dphi_lepcuts[-2],None)}
        #########################################
        #outdir = "./combine/HNL_datacards/muon_v11/"   ### v11 , full run 2, looseeID, new timing 
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v12/"   ### v12 , full run 2,no looseID, new timing 
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v13/"   ### v13 , full run 2,no looseID, new timing, no mu pT cut for signal
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v14/"   ### v14 , new xSec, new trig. bits 
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v15/"   ### v15 , fixed muon SF, OOTpred 
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v16/"   ### v16 , as in v15 but , high nHit T.F. for DT
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v17/"   ### v17 , as in v15, 150 nHit for DT(Dec23 signal), OOTxTF 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v18/"   ### v18 , as in v17, (Dec23 signal), ABCD 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #outdir = "./combine/HNL_datacards/muon_v19/"   ### v19 , as in v17, (Jan9 signal), ABCD, split MB2,MB34, apply CSC,DT_cls_SF 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        outdir = "./combine/HNL_datacards/muon_v26_full_lumi/"   ### v20 , as in v19, (Jan9 signal), unblind 
        cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #########################################
        #outdir = "./combine/HNL_datacards/tau_v1/mu/"   ### v1 
        #cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        #outdir = "./combine/HNL_datacards/tau_v2/mu/"    ### v2 , use real ABC, 150 DT, more ctau points 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #outdir = "./combine/HNL_datacards/tau_v3/mu/"    ### v3 (as in v18), split MB2,MB34, apply CSC,DT cls SF 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        #outdir = "./combine/HNL_datacards/tau_v4/mu/"    ### v4 (as in v20), unblind 
        #cut = {"CSC":(200,2.8,None), "DT":(150,2.8,None)}
        isMuon=True

        f_yield = outdir+"yields.json"
        
        print("Working on muon Channel: ")
        print("        yields = ", f_yield)
        if options.muon:
            if options.writeYields: writeYields(cut,isMuon,f_yield,True,shifts,options.tauSignals,options.unblind,options.useOOT) 
            else:
                makeAllcards(f_yield.replace(".json","_mixed.json"),outdir,"",options.dryRun)
                #makeAllcards(f_yield.replace(".json","_dirac.json"),outdir,"",options.dryRun)

