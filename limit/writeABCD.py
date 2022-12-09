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
    text_file.write('------------------------------ \n')

def writeUncSection(text_file,signal_rate,bkg_proc,lnN_unc):
    ### unc = {"proc_uncX":[a,b,c,d]}

    ## match num. of proc
    procs = {} # name:[a,b,c,d]
    procs.update(signal_rate)
    procs.update(bkg_proc)

    for unc_name,unc_arr in lnN_unc.items():
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
    text_file.write('norm  rateParam       *     ggH      1 \n')
    #### uncertainties ####
    writeUncSection(text_file,signal_rate,bkg_rate,bkg_unc)


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
def calSigal(muon=False):
    if muon:
        data = {
        # from ct1m ###"HNL_muonType_mHNL2p0_pl15000": {'CSC': np.array([ 1.23607486,  2.94705472, 33.51721181, 10.69697318]), 'DT': np.array([ 2.34226239,  2.01891048, 30.98168154, 19.26501103]), 'norm': 1},
        ##from ct1m ###"HNL_muonType_mHNL2p0_pl2000": {'CSC': np.array([ 0.7164319 ,  1.72425638, 19.74654625,  6.25635821]), 'DT': np.array([ 1.37572843,  1.19858567, 18.42712308, 11.41849004]), 'norm': 1},
        #"HNL_muonType_mHNL1p0_pl10" :{'CSC': np.array([1313.72597472, 2083.29270286, 2387.01682551, 1588.47727032]), 'DT': np.array([1245.3276547 ,    0.        ,  514.90494597,    0.        ]), 'norm': 100},
        #"HNL_muonType_mHNL1p0_pl100" :{'CSC': np.array([  958.36080102,  2541.21051879, 35942.66029522, 15694.03822141]), 'DT': np.array([ 2516.41081957,  2687.36740861, 24510.1580307 , 20040.76760296]), 'norm': 100},
        #"HNL_muonType_mHNL1p0_pl1000": {'CSC': np.array([  28.82915809,   72.20854441, 1063.38960708,  482.29479901]), 'DT': np.array([ 55.81080617,  84.18083673, 980.13814336, 767.72333672]), 'norm': 10},
        #"HNL_muonType_mHNL1p0_pl2000": {'CSC':np.array([  7.5398423 ,  19.03590541, 286.4139182 , 127.38951007]), 'DT': np.array([ 14.8699539 ,  22.64145349, 266.42989897, 207.79041764]), 'norm': 10},
        #"HNL_muonType_mHNL1p0_pl5000": {'CSC':np.array([ 1.23971623,  3.1465619 , 48.01300414, 21.06873556]), 'DT': np.array([ 2.47262373,  3.78655542, 44.83523666, 34.87615238]), 'norm': 1},
        #"HNL_muonType_mHNL1p0_pl10000": {'CSC': np.array([ 0.        ,  2.25792804, 12.5938797 ,  3.91385083]), 'DT': np.array([ 0.54224782,  1.20117507, 13.53483871,  7.56249241]), 'norm': 1},
        #"HNL_muonType_mHNL1p0_pl20000": {'CSC': np.array([0.        , 0.58674456, 3.23417357, 0.98425379]), 'DT': np.array([0.13650751, 0.3022155 , 3.412853  , 1.90518924]), 'norm': 1},
        #"HNL_muonType_mHNL1p0_pl50000": {'CSC': np.array([0.        , 0.0962575 , 0.52797146, 0.15803935]), 'DT': np.array([0.02193248, 0.04854013, 0.54887452, 0.30623774]), 'norm': 1},
        #"HNL_muonType_mHNL2p0_pl10": {'CSC': np.array([12.94954136,  0.        ,  8.14461429,  0.        ]), 'DT': np.array([ 2.89209086,  0.        , 12.84317809,  0.        ]), 'norm': 1},
        #"HNL_muonType_mHNL2p0_pl100": {'CSC': np.array([ 46.08049672,  65.31234612, 763.99850184, 382.83887428]), 'DT': np.array([ 66.71386847,  59.28390516, 373.94752874, 309.02756679]), 'norm': 1},
        #"HNL_muonType_mHNL2p0_pl1000": {'CSC': np.array([ 2.62177929,  6.13920546, 68.80480406, 22.27152984]), 'DT': np.array([ 4.84003038,  4.08741815, 62.3940202 , 39.06735284]), 'norm': 1},
        #"HNL_muonType_mHNL2p0_pl10000":{'CSC': np.array([0.01535685, 0.08224103, 1.12456859, 0.30455309]), 'DT': np.array([0.06101222, 0.04465391, 0.90300775, 0.63003223]), 'norm': 1},
        #"HNL_muonType_mHNL2p0_pl2000" :{'CSC': np.array([ 0.35527838,  1.80945444, 25.12737918,  6.99742389]), 'DT': np.array([ 1.35857254,  1.04874446, 19.70587762, 13.900502  ]), 'norm': 1},
        #"HNL_muonType_mHNL2p0_pl5000" :{'CSC': np.array([0.06024815, 0.31859128, 4.37329474, 1.19273078]), 'DT': np.array([0.23706951, 0.17584326, 3.49085573, 2.44230647]), 'norm': 1},
        #"HNL_muonType_mHNL2p0_pl800": {'CSC': np.array([  3.92207907,   9.06671745, 100.43177256,  32.844178  ]), 'DT': np.array([ 7.10264082,  5.91173037, 89.76806317, 56.49657938]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl50": {'CSC': np.array([0.52972903, 0.41909375, 1.16959573, 1.26086066]), 'DT': np.array([0.38047228, 0.29476375, 0.12242464, 0.05909414]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl60": {'CSC': np.array([0.68995215, 0.53018717, 2.16536157, 2.10749023]), 'DT': np.array([0.48909917, 0.42485259, 0.30484776, 0.16079799]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl80": {'CSC': np.array([0.9071468 , 0.73450322, 4.61030034, 3.84690764]), 'DT': np.array([0.71360755, 0.67620385, 0.99318524, 0.5558966 ]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl90": {'CSC': np.array([0.96786351, 0.81570234, 5.85914099, 4.60951678]), 'DT': np.array([0.81912794, 0.78479803, 1.46442202, 0.82713588]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl100": {'CSC': np.array([1.00382177, 0.88081391, 7.03360537, 5.26295055]), 'DT': np.array([0.91579247, 0.87776941, 1.97958913, 1.12285429]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl200": {'CSC': np.array([ 0.89609229,  1.07542741, 13.89661744,  6.37708174]), 'DT': np.array([0.84830022, 1.07858285, 7.20642776, 4.46332663]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl500": {'CSC': np.array([0.32837435, 0.56626363, 9.24723856, 3.59179129]), 'DT': np.array([0.43360458, 0.75001888, 7.07675428, 3.89266032]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl800": {'CSC': np.array([0.16322816, 0.31673343, 5.34028832, 1.96939856]), 'DT': np.array([0.23616548, 0.44824362, 4.48827116, 2.38402303]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl1000": {'CSC': np.array([0.11367208, 0.23031757, 3.90969984, 1.4149658 ]), 'DT': np.array([0.16994022, 0.33335487, 3.38758482, 1.77729034]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl1500": {'CSC': np.array([0.05674916, 0.12227352, 2.08696206, 0.7353985 ]), 'DT': np.array([0.08887886, 0.18265597, 1.8813124 , 0.97029353]), 'norm': 1},
        #"HNL_muonType_mHNL4p0_pl2000": {'CSC': np.array([0.03388892, 0.07544135, 1.28896805, 0.44779054]), 'DT': np.array([0.05441054, 0.11464907, 1.18446408, 0.60551109]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl5000":{'CSC': np.array([0.00073515, 0.00511466, 0.14058247, 0.04399492]), 'DT': np.array([0.00582457, 0.01495561, 0.14932642, 0.06974714]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl1000": {'CSC': np.array([0.05042282, 0.13644857, 2.27044505, 0.85899962]), 'DT': np.array([0.12006929, 0.2090709 , 1.91663381, 1.03892495]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl100" :{'CSC': np.array([0.33456926, 0.54080644, 3.02210338, 1.83951148]), 'DT': np.array([0.52526878, 0.29551597, 0.66847637, 0.3938327 ]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl200" :{'CSC': np.array([0.33693308, 0.63268827, 6.82774689, 3.33129667]), 'DT': np.array([0.59264865, 0.57457409, 3.13038029, 1.89648114]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl50"  :{'CSC': np.array([0.12003204, 0.25953659, 0.33222411, 0.2711571 ]), 'DT': np.array([0.24203619, 0.05143207, 0.02458772, 0.00793734]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl60" :{'CSC': np.array([0.18028231, 0.33683376, 0.7032576 , 0.5298672 ]), 'DT': np.array([0.32015663, 0.09241558, 0.06911001, 0.02999097]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl70" :{'CSC': np.array([0.23361351, 0.40223641, 1.19865222, 0.84395739]), 'DT': np.array([0.38731178, 0.14067563, 0.15231426, 0.0768085 ]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl200" : {'CSC': np.array([0.33693308, 0.63268827, 6.82774689, 3.33129667]), 'DT': np.array([0.59264865, 0.57457409, 3.13038029, 1.89648114]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl500" : {'CSC': np.array([0.14162324, 0.33918348, 5.14047411, 2.09490877]), 'DT': np.array([0.30551087, 0.44207803, 3.70918989, 2.10335137]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl800" : {'CSC': np.array([0.07199808, 0.18838995, 3.0681604 , 1.18362337]), 'DT': np.array([0.16697468, 0.27639257, 2.48869457, 1.36663156]), 'norm': 1},
        #"HNL_muonType_mHNL4p5_pl5000": {'CSC': np.array([0.00073515, 0.00511466, 0.14058247, 0.04399492]), 'DT': np.array([0.00582457, 0.01495561, 0.14932642, 0.06974714]), 'norm': 1},
        #"HNL_muonType_mHNL5p0_pl100": {'CSC': np.array([0.34527562, 0.32817664, 1.06063231, 0.70045599]), 'DT': np.array([0.15490411, 0.10799062, 0.28216942, 0.17626673]), 'norm': 1},
        #"HNL_muonType_mHNL5p0_pl50" : {'CSC': np.array([0.1857032 , 0.11720223, 0.10026844, 0.0758087 ]), 'DT': np.array([0.07068723, 0.01439644, 0.08330753, 0.00342477]), 'norm': 1},
        #"HNL_muonType_mHNL5p0_pl200": {'CSC': np.array([0.28719842, 0.38748771, 2.99977646, 1.52883504]), 'DT': np.array([0.2140579 , 0.23215314, 1.42081393, 0.98332443]), 'norm': 1},
        #"HNL_muonType_mHNL5p0_pl500": {'CSC': np.array([0.11016925, 0.21612216, 2.6892949 , 1.10706214]), 'DT': np.array([0.14367971, 0.19494229, 2.07599607, 1.25006474]), 'norm': 1},
        #"HNL_muonType_mHNL5p0_pl800": {'CSC': np.array([0.05533046, 0.12587289, 1.68748153, 0.65352   ]), 'DT': np.array([0.08575991, 0.12676015, 1.48230975, 0.84895672]), 'norm': 1},
        #"HNL_muonType_mHNL5p0_pl1000": {'CSC': np.array([0.03866598, 0.0936598 , 1.27124608, 0.48186906]), 'DT': np.array([0.06364373, 0.09759387, 1.16621215, 0.65594387]), 'norm': 1},
        #"HNL_muonType_mHNL5p0_pl5000": {'CSC': np.array([0.00184826, 0.00352256, 0.09103392, 0.0386287 ]), 'DT': np.array([0.00481397, 0.00941965, 0.10353558, 0.04951436]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl115" :{'CSC': np.array([0.42879209, 0.3762481 , 3.0044719 , 2.24811973]), 'DT': np.array([0.39118953, 0.37494761, 0.84560046, 0.47963797]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl118" :{'CSC': np.array([0.377185  , 0.33096492, 2.64286999, 1.97754826]), 'DT': np.array([0.34410808, 0.32982095, 0.74382858, 0.42191135]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl1150":{'CSC': np.array([0.04855612, 0.09838235, 1.67006573, 0.6044162 ]), 'DT': np.array([0.07259159, 0.14239572, 1.4470393 , 0.75918659]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl1175":{'CSC': np.array([0.04271217, 0.08654159, 1.4690657 , 0.53167195]), 'DT': np.array([0.06385486, 0.12525775, 1.27288152, 0.66781502]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl1725":{'CSC': np.array([0.02424095, 0.05223031, 0.89146583, 0.31413251]), 'DT': np.array([0.03796545, 0.07802325, 0.80362061, 0.41447018]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl1762":{'CSC': np.array([0.02132344, 0.04594415, 0.78417385, 0.27632523]), 'DT': np.array([0.03339614, 0.06863279, 0.7069012 , 0.3645868 ]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl230" :{'CSC': np.array([0.3827744 , 0.45937912, 5.93607324, 2.72403154]), 'DT': np.array([0.36235956, 0.46072699, 3.07829464, 1.90655272]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl235" :{'CSC': np.array([0.33670576, 0.40409074, 5.22163972, 2.39618191]), 'DT': np.array([0.31874794, 0.4052764 , 2.70780783, 1.6770904 ]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl2300":{'CSC': np.array([0.01447597, 0.0322255 , 0.55059505, 0.19127802]), 'DT': np.array([0.02324198, 0.04897345, 0.50595518, 0.25864986]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl2350":{'CSC': np.array([0.01273373, 0.02834701, 0.48432842, 0.16825684]), 'DT': np.array([0.0204447 , 0.04307927, 0.44506116, 0.22752017]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl575" :{'CSC': np.array([0.14026825, 0.24188493, 3.95004652, 1.53426805]), 'DT': np.array([0.18521835, 0.32037775, 3.02290337, 1.66278714]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl588" :{'CSC': np.array([0.12338633, 0.21277298, 3.47464039, 1.34961189]), 'DT': np.array([0.16292648, 0.28181882, 2.65908315, 1.46266312]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl920" :{'CSC': np.array([0.06972447, 0.13529572, 2.2811553 , 0.84124746]), 'DT': np.array([0.10088035, 0.19147155, 1.91720801, 1.01835827]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl940" :{'CSC': np.array([0.06133281, 0.11901226, 2.00660785, 0.73999949]), 'DT': np.array([0.08873894, 0.16842708, 1.68646328, 0.8957942 ]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl57"  :{'CSC': np.array([0.22627883, 0.17901991, 0.49960402, 0.53858871]), 'DT': np.array([0.16252238, 0.12591116, 0.05229486, 0.02524263]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl59"  :{'CSC': np.array([0.19904514, 0.15747405, 0.43947439, 0.4737671 ]), 'DT': np.array([0.14296207, 0.11075718, 0.04600093, 0.02220456]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl92"  :{'CSC': np.array([0.38749644, 0.31375008, 1.96933395, 1.64324344]), 'DT': np.array([0.3048243 , 0.28884695, 0.42424859, 0.23745655]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl94"  :{'CSC': np.array([0.34085948, 0.27598883, 1.73231562, 1.44547159]), 'DT': np.array([0.2681373 , 0.2540829 , 0.37318833, 0.20887757]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl69"  :{'CSC': np.array([0.29471967, 0.22647453, 0.92495493, 0.90023464]), 'DT': np.array([0.20892339, 0.18147985, 0.13021864, 0.0686864 ]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl70"  :{'CSC': np.array([0.25924881, 0.19921728, 0.81363239, 0.79188729]), 'DT': np.array([0.18377851, 0.15963792, 0.11454623, 0.06041968]), 'norm': 1},
        "HNL_muonType_mHNL4p6_pl103" :{'CSC': np.array([0.41343217, 0.34843507, 2.50278819, 1.96899924]), 'DT': np.array([0.34989835, 0.335234  , 0.62554189, 0.35331901]), 'norm': 1},
        "HNL_muonType_mHNL4p7_pl106" :{'CSC': np.array([0.36367372, 0.30649931, 2.20156621, 1.73202119]), 'DT': np.array([0.30778648, 0.29488705, 0.55025507, 0.31079546]), 'norm': 1},
        }
    else:
        data={
        "HNL_electronType_mHNL10p0_pl100"  :{'CSC': np.array([0.00234753, 0.00051788, 0.        , 0.        ]), 'DT': np.array([0.0011039, 0.       , 0.       , 0.       ]), 'norm': 1}                     ,
        "HNL_electronType_mHNL10p0_pl1000" :{'CSC': np.array([0.00391218, 0.00314663, 0.00524341, 0.00604237]), 'DT': np.array([0.00192836, 0.0017127 , 0.00892355, 0.00903163]), 'norm': 1}                 ,
        "HNL_electronType_mHNL10p0_pl10000":{'CSC': np.array([0.00011476, 0.00013608, 0.00023872, 0.0001808 ]), 'DT': np.array([4.06912175e-05, 5.08776222e-05, 4.30198853e-04, 3.17931076e-04]), 'norm': 1} ,
        "HNL_electronType_mHNL1p0_pl1000"  :{'CSC': np.array([ 84.32759251, 111.64981372, 156.62843274, 168.77611882]), 'DT': np.array([ 46.30755564,  48.15297065, 322.25762851, 369.49998721]), 'norm': 1} ,
        "HNL_electronType_mHNL1p0_pl10000" :{'CSC': np.array([0.50949814, 3.36698952, 1.59614641, 2.25987333]), 'DT': np.array([0.66427956, 0.        , 3.64510972, 4.06372616]), 'norm': 1}                 ,
        "HNL_electronType_mHNL1p0_pl20000" :{'CSC': np.array([0.12863117, 0.92021893, 0.40821676, 0.56840318]), 'DT': np.array([0.16727466, 0.        , 0.91957047, 1.02408306]), 'norm': 1}                 ,
        "HNL_electronType_mHNL1p0_pl50000" :{'CSC': np.array([0.02070258, 0.16096901, 0.06627455, 0.09127591]), 'DT': np.array([0.02688028, 0.        , 0.14793334, 0.16464103]), 'norm': 1}                 ,
        "HNL_electronType_mHNL2p0_pl10000" :{'CSC': np.array([0.08077014, 0.0511257 , 0.2136278 , 0.12036926]), 'DT': np.array([0.02120186, 0.04912579, 0.17766103, 0.19017568]), 'norm': 1}                 ,
        "HNL_electronType_mHNL2p0_pl1000"  :{'CSC': np.array([ 5.8418507 ,  4.10399175, 15.52154069,  9.89584296]), 'DT': np.array([ 1.60758728,  2.94298224, 12.68235199, 14.58692474]), 'norm': 1}         ,
        "HNL_electronType_mHNL2p0_pl2000"  :{'CSC': np.array([1.74147344, 1.15915676, 4.44141076, 2.75812382]), 'DT': np.array([0.46869685, 0.95166955, 3.81258008, 4.22303001]), 'norm': 1}                 ,
        "HNL_electronType_mHNL2p0_pl5000"  :{'CSC': np.array([0.31115517, 0.19956601, 0.79587618, 0.47109629]), 'DT': np.array([0.08223913, 0.18358378, 0.6837058 , 0.73842168]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl100"   :{'CSC': np.array([0.94907259, 0.67067039, 0.86846901, 1.17488095]), 'DT': np.array([0.59362673, 0.44108409, 0.54051344, 0.29325154]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl1000"  :{'CSC': np.array([0.28226471, 0.28934621, 0.50104768, 0.48555935]), 'DT': np.array([0.1330846 , 0.15056762, 1.02402093, 0.91202053]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl200"   :{'CSC': np.array([1.59274261, 0.94461057, 1.63835251, 1.92739285]), 'DT': np.array([0.71717849, 0.65136852, 2.18898571, 2.05707483]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl2000"  :{'CSC': np.array([0.08730201, 0.0990024 , 0.1663454 , 0.15538752]), 'DT': np.array([0.04320927, 0.04948626, 0.35750548, 0.3149966 ]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl500"   :{'CSC': np.array([0.75124284, 0.65635638, 1.16506758, 1.19958757]), 'DT': np.array([0.33695848, 0.36669251, 2.14576641, 1.94433359]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl800"   :{'CSC': np.array([0.39742771, 0.38986524, 0.68165999, 0.67159408]), 'DT': np.array([0.18413363, 0.20666247, 1.3578237 , 1.21508386]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl50"    :{'CSC': np.array([0.35465326, 0.2493752 , 0.09292032, 0.15007424]), 'DT': np.array([0.1730583 , 0.12490802, 0.02713038, 0.01072008]), 'norm': 1}                 ,
        "HNL_electronType_mHNL4p0_pl80"    :{'CSC': np.array([0.72826691, 0.50013602, 0.50762429, 0.72367099]), 'DT': np.array([0.43639043, 0.32789811, 0.26002246, 0.1317294 ]), 'norm': 1}                 ,
        "HNL_electronType_mHNL7p0_pl100"   :{'CSC': np.array([0.02695604, 0.01149086, 0.00645032, 0.0145628 ]), 'DT': np.array([0.01397228, 0.00687955, 0.        , 0.00148831]), 'norm': 1}                 ,
        "HNL_electronType_mHNL7p0_pl1000"  :{'CSC': np.array([0.02228937, 0.02672522, 0.0379208 , 0.04380189]), 'DT': np.array([0.01006018, 0.00889886, 0.06141662, 0.06364733]), 'norm': 1}                 ,
        "HNL_electronType_mHNL7p0_pl10000" :{'CSC': np.array([0.00054382, 0.00054071, 0.00085457, 0.00081199]), 'DT': np.array([2.99369415e-04, 9.73379035e-05, 2.03206278e-03, 1.87740821e-03]), 'norm': 1} ,
        }
    for d in data:
        print(d,data[d])
    return data

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

## hard-coded scale for muon channel
def scale2017(h1,h2):
    ##2016+2018 = 51.2+30.9  =82.1## mu
    h1.scale(82.1)
    ##2017 = 37.9 ## mu
    h2.scale(37.9)
    h1.add(h2)      ## this is done in-place
    return 

def writeYields(cut = None,muon=True,outf="yields.json",debug=True,shifts=None):
    from coffea import hist
    if cut==None:
        ## nHit, dphi_lep, dphi_MET
        dphi_lepcuts = np.linspace(0,np.pi,31)[1:-2]
        if muon:
            cut = {"CSC":(220,dphi_lepcuts[-10],0.7), "DT":(110,dphi_lepcuts[-8],0.7)}
        else:
            cut = {"CSC":(160,dphi_lepcuts[-2],0.7) , "DT":(100,dphi_lepcuts[-10],0.7)}

    if muon:
        lumi = 120
        #bkg     = loadbkg("../HNL_histograms_Feb23_muons_data.pickle",True,cut,debug)
        #signals = loadhist("../HNL_histograms_Feb23_muons_signal.pickle",True,cut,debug)
        #signals.update( loadhist("../HNL_histograms_Mar1_muons_signal.pickle",True,cut,debug))
        #bkg     = loadbkg("../HNL_histograms_Nov29_muons_all.pickle",True,cut,debug)
        bkg      = loadbkg("../HNL_histograms_Dec6_muons_data.pickle",True,cut,debug)
        signals  = loadhist("../HNL_histograms_Dec6_muons_signals.pickle",True,cut,debug,lumi,
                            "../HNL_histograms_Dec6_muons_signals_2017.pickle" )
        #signals = loadSignalFromJson("./limitpoints_muon.json",True,cut,debug,lumi) ## this is too slow 
    else:
        lumi = 122
        #bkg     = loadbkg("../HNL_histograms_Feb3_electrons.pickle",False,cut,debug)
        #signals = loadhist("../HNL_histograms_Apr1_ele_signal.pickle",False,cut,debug)
        #bkg     = loadbkg("../HNL_histograms_Nov29_ele_all.pickle",False,cut,debug)
        bkg     = loadbkg("../HNL_histograms_Dec6_ele_data.pickle",False,cut,debug)
        signals = loadhist("../HNL_histograms_Dec6_ele_signals.pickle",False,cut,debug,lumi)
        #signals = loadSignalFromJson("./limitpoints_ele.json",False,cut,debug,lumi)

    #data = {**bkg,**signals} 
    data = {}
    for k,v in bkg.items(): data[k] = v 
    for k,v in signals.items(): data[k] = v 
    for k,v in data.items():
       v['DT'] = v['DT'].tolist()
       v['CSC'] = v['CSC'].tolist()
       v['CSC_unc'] = v['CSC_unc'].tolist()
       v['DT_unc'] = v['DT_unc'].tolist()
    with open(outf,"w") as f:
       f.write(json.dumps(data,indent=4))
    if shifts is not None:
        shiftYields(outf,shifts)
    return 

def loadbkg(fin='../HNL_histograms_Feb3_electrons.pickle',muon=False,cut=None,debug=True):
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
    OOT = h.integrate("region","ABCD_OOT")
    hdt = out['dphi_cluster_dt'].integrate("dataset",datasets)
    OOT_dt = hdt.integrate("region","ABCD_dt_OOT")

    if cut==None:
        ## nHit, dphi_lep, dphi_MET
        dphi_lepCuts = np.linspace(0,np.pi,31)[1:-2]
        # dphi_lepCuts[-2] ### 2.827
        # dphi_lepCuts[-10] ###. 1.9896
        if muon:
            cut = {"CSC":(160,dphi_lepCuts[-2],0.7), "DT":(100,1.989,0.7)}
        else:
            cut = {"CSC":(160,dphi_lepCuts[-2],0.7), "DT":(100,dphi_lepCuts[-10],0.7)}
    kfactor=0.25
    CSC,CSC_unc = predIntimeFromOOT(OOT,cut['CSC'][0],cut["CSC"][1],cut['CSC'][2],False,kfactor,lumiScale)
    kfactor=0.9
    DT,DT_unc = predIntimeFromOOT(OOT_dt,cut["DT"][0],cut["DT"][1],cut["DT"][2],False,kfactor,lumiScale)

    if debug:
        print("bkg CSC = " ,CSC)
        print("bkg CSCunc = " ,CSC_unc)
        print("bkg DT = " ,DT)
        print("bkg DTunc = ", DT_unc)
    data={}
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
                print(f"     processing time =   {toc - tic:0.4f} seconds")
            del(events)
            data[dataset] = scaleHist(out,dataset,lumi,cut,debug) 
    stop = time.perf_counter()
    if debug: print(f"     total processing time =   {stop - start:0.4f} seconds")
    return data

    
def loadhist(fin='../HNL_histograms_Feb23_muons_signal.pickle',muon=True,cut=None,debug=True,lumi=137,fin2017=None):
    from coffea import hist
    
    with open(fin,'rb')  as f:
        out = pickle.load(f)
   
    
    if fin2017 is not None:
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

        dphi_lepCuts = np.linspace(0,np.pi,31)[1:-2]
        ### NOT USED FOR SIGNAL
        lumi = 1
        kfactor=1     ## Muon, CSC InT/OOT 

        ## CSC cuts
        CSC,CSC_unc = predIntimeFromOOT(signal,cut["CSC"][0],cut["CSC"][1],cut["CSC"][2],True,kfactor,lumi)

        ### DT cuts
        DT,DT_unc = predIntimeFromOOT(signal_dt,cut["DT"][0],cut["DT"][1],cut["DT"][2],True,kfactor,lumi)
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
        #data[sample_name] = {"CSC":CSC,"DT":DT,"norm":1}
        data[sample_name] = {"CSC":CSC,"DT":DT,"CSC_unc":CSC_unc,"DT_unc":DT_unc,"norm":1}
              
    return data 

def makeAllcards(f_yield,outdir="./combine/HNL_datacards/",suffix="",test=False):

    if not os.path.exists(outdir):
        print("mk dir = ", outdir)
        os.mkdir(outdir)
    #outdir = "../combine/dt_datacards/"
    #outdir = "./combine/HNL_datacards/"
    norm = 100.0    # 1% BR
    suffix = ""


    #ZmumuCR = CSC:73      DT:172  (200,130)
    #ZmumuCR = CSC:54      DT:130  (220,150)
    

    bkg_proc_CSC= {
     #   "Zmumu_CSC":[-1,-1,-1,2.07], ## 2.5% TF 
        #"Zmumu_CSC":[-1,-1,-1,4.68], ## 5.6% TF 
        "Zmumu_CSC":[-1,-1,-1,3.02], ## 5.6% TF, 220
    }
    bkg_proc_DT= {
     #   "Zmumu_DT":[-1,-1,-1,1.95],    # 1.0% TF
        #"Zmumu_DT":[-1,-1,-1,2.56],     # 1.3% TF
        "Zmumu_DT":[-1,-1,-1,1.69],     # 1.3% TF, 150
    }
    # procName_uncName
    bkg_unc = {
        "bkg_sys_D":[0  ,0 ,0 ,0.5],
        #"Zmumu_CSC_sys_D":[0 ,0 ,0 ,0.32],      
        #"Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.2 ],
        "Zmumu_CSC_sys_D":[0 ,0 ,0 ,0.25],    
        "Zmumu_DT_sys_D" :[0 ,0 ,0 ,0.23],
    }
    sig_unc = {
        "clus":{
            "ggH":[0.15 ,0.15, 0.15 ,0.15],
        },
        "lumi":{
            "ggH":[0.018 ,0.018, 0.018 ,0.018],
        },
    }
    with open(f_yield,'r') as f:
        data = json.load(f)
        bkg_rate_CSC = {}       # {bkg: [a,b,c,d]}
        bkg_rate_CSC.update({"bkg":data['bkg']["CSC"]})          # {bkg: [a,b,c,d]}
        bkg_rate_CSC.update(bkg_proc_CSC)                        # add other bkg 
        bkg_rate_DT = {}       # {bkg: [a,b,c,d]}
        bkg_rate_DT.update({"bkg":data['bkg']["DT"]})           # {bkg: [a,b,c,d]}
        bkg_rate_DT.update(bkg_proc_DT)                         # add other bkg 
        for name,signal in data.items():
            #if not ("3p1" in name or "3p2" in name):  continue
            if suffix:
                name = name+suffix
            #norm = 1
            norm = signal["norm"]
            sigRate = {"ggH":np.array(signal["CSC"])/norm }
            obs = bkg_rate_CSC['bkg']
            obs[-1] = bkg_rate_CSC['bkg'][0]*bkg_rate_CSC['bkg'][2]/bkg_rate_CSC['bkg'][1]  ## force D=A*C/B
            make_datacard_2sig(outdir,name+"_CSC", sigRate, norm, bkg_rate_CSC, obs, bkg_unc,  sig_unc)
            sigRate = {"ggH":np.array(signal["DT"]) /norm}
            obs = bkg_rate_DT['bkg']
            obs[-1] = bkg_rate_DT['bkg'][0]*bkg_rate_DT['bkg'][2]/bkg_rate_DT['bkg'][1]  ## force D=A*C/B
            make_datacard_2sig(outdir,name+"_DT", sigRate, norm, bkg_rate_DT, obs, bkg_unc,  sig_unc)
   
            def Run(cmd,test=False):
                print(cmd)
                if not test: os.system(cmd)
            csc_limit = "combine -M AsymptoticLimits {odir}{name}_CSC.txt -n _{name}_CSC --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
            
            Run(csc_limit,test)
            Run("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_CSC",outdir),test)

            dt_limit = "combine -M AsymptoticLimits {odir}{name}_DT.txt -n _{name}_DT --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
            Run(dt_limit,test)
            Run("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_DT",outdir),test)

            cmd = "python combination.py CSC={odir}{name}_CSC.txt DT={odir}{name}_DT.txt {odir}{name}_comb.txt".format(name=name,odir=outdir)
            Run(cmd,test)
            cmd = "combine -M AsymptoticLimits {odir}{name}_comb.txt -n _{name}_comb --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
            Run(cmd,test)
            Run("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_comb",outdir),test)
               
import sys
sys.path.insert(0,"../")
#from HNLprocessor.util import f_1m   ## Move to util if possible
def f_1m(x):
    x0 = np.array([1,2,4,5,7,10])
    y0 = np.array([13.57,0.4238,0.01289,0.004156,0.0007452,0.000121])    
    return np.exp(np.poly1d(np.polyfit(x0,np.log(y0),5))(x))

def f_xsec(m):
    def xsec_m(x):
        return f_1m(m)/(x/1000.)
    return xsec_m 

# Find the corresponding yield of (m,ct)-> (m,xx)
def shift_ctau(N_yield,m_old,ct_old,m_new):
    ct_new = (m_new/m_old) * ct_old
    return N_yield * f_xsec(m_new)(ct_new)/f_xsec(m_old)(ct_old)

# shift yielsd in f_yield.json according to shifts[{"m_src":int,"m_target":int}]
def shiftYields(f_yield,shifts):
    ## init
    for shift in shifts:
        shift['CSC_target']=[]
        shift['DT_target']=[]
        shift['CSC_unc_target']=[]
        shift['DT_unc_target']=[]
        shift['names']=[]
    with open(f_yield,'r') as f:
        data = json.load(f)
        for name,signal in data.items():
            if name=="bkg": continue
            m =float(name.split("_")[-2].replace("mHNL","").replace("p","."))
            ct =float(name.split("_")[-1].replace("pl",""))    
            for shift in shifts:
                if m==shift["m_src"]:     
                    ct_new = (shift['m_target']/m) * ct            
                    shift["CSC_target"].append(shift_ctau(np.array(signal['CSC']),m,ct,shift['m_target']))            
                    shift["DT_target"].append(shift_ctau(np.array(signal['DT']),m,ct,shift['m_target']))                        
                    shift["CSC_unc_target"].append(shift_ctau(np.array(signal['CSC_unc']),m,ct,shift['m_target']))            
                    shift["DT_unc_target"].append(shift_ctau(np.array(signal['DT_unc']),m,ct,shift['m_target']))                        
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
                    "DT":shift["DT_target"][i].tolist(),  
                    "CSC_unc":shift["CSC_unc_target"][i].tolist(),
                    "DT_unc":shift["DT_unc_target"][i].tolist(),  
                    "norm":1,
                }  
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
    parser.add_option('--muon', dest='muon', action='store_true',default = False, help='make muon datacard')
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
        cut = {"CSC":(200,dphi_lepcuts[-4],None), "DT":(130,dphi_lepcuts[-4],None)}
        isMuon=True
        #f_yield = "yields.json"
        f_yield = "./combine/HNL_datacards/muon_v6/yields.json"
        makeAllcards(f_yield,outdir,"",True)
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
        outdir = "./combine/HNL_datacards/ele_v8/"   ###  v8, full run 2, looseID, new timing, 1 muon pT cut
        cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        isMuon=False

        f_yield = outdir+"yields.json"
        shifts = [
            {"m_src":4.0,"m_target":3.0},
            {"m_src":4.0,"m_target":3.1},
            {"m_src":4.0,"m_target":3.2},
            {"m_src":4.0,"m_target":3.5},
        ]
        if not options.muon:
            if options.writeYields: writeYields(cut,isMuon,f_yield,True,shifts) 
            else:   makeAllcards(f_yield,outdir,"",options.dryRun)
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
        print("Working on muon Channel: ")
        outdir = "./combine/HNL_datacards/muon_v11/"   ### v11 , full run 2, loosID, new timing 
        cut = {"CSC":(200,2.8,None), "DT":(130,2.8,None)}
        isMuon=True

        f_yield = outdir+"yields.json"
        if options.muon:
            if options.writeYields: writeYields(cut,isMuon,f_yield,True,shifts) 
            else:            makeAllcards(f_yield,outdir,"",options.dryRun)
