import glob,os
import json 

fileset = {
     #r'HNL,$m_N$=5': [
     #    '/eos/uscms/store/user/kkwok/llp/HNL/HeavyNeutralLepton_Tree.root',
     #],

#         'EGamma_2018A': glob.glob("root://cmsxrootd.fnal.gov//store/lpclonglived/HNL/EGamma_2018A/HeavyNeutralLepton_Tree_*.root"),
    
      'EGamma_2018A': glob.glob("/eos/uscms/store/user/lpclonglived/HNL/EGamma_2018A/HeavyNeutralLepton_Tree_*.root"),
    
#      'EGamma_2018A': glob.glob("/uscms/home/kkwok/lpclonglived/HNL/EGamma_2018A/HeavyNeutralLepton_Tree_*.root"),
#       'EGamma_2018A': ["/uscms/home/kkwok/lpclonglived/HNL/EGamma_2018A/HeavyNeutralLepton_Tree_0.root"],
    
#     r'$m_s$=15,$c\tau$=1m': [
#         "~/eos/llp/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/HeavyNeutralLepton_Tree_15_1000.root"
#     ],    
#     'WJetsToLNu': [
#          "~/eos/llp/WJetsToLNu/HeavyNeutralLepton_Tree_0.root"        
#     ],    
       'WJetsToLNu':glob.glob("/eos/uscms/store/user/kkwok/llp/WJetsToLNu/HeavyNeutralLepton_Tree_*.root")        
}

signals = {
     #'HNL_testpoint1': [
     #    '/eos/uscms/store/user/kkwok/llp/HNL/HeavyNeutralLepton_Tree.root',
     #],
     "HNL_electronType_mHNL1p0_pl10"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL1p0_pl10/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL1p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL1p0_pl100/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL1p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL1p0_pl1000/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL2p0_pl10"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL2p0_pl10/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL2p0_pl10_rwctau4"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL2p0_pl10/HeavyNeutralLepton_Tree.root"),
     ##"HNL_electronType_mHNL2p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL2p0_pl100/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL2p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL2p0_pl1000/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL2p0_pl1000_rwctau10000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL2p0_pl1000/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL4p0_pl10"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL4p0_pl10/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL4p0_pl10_rwctau40"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL4p0_pl10/HeavyNeutralLepton_Tree.root"),
     ##"HNL_electronType_mHNL4p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL4p0_pl100/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL4p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL4p0_pl1000/HeavyNeutralLepton_Tree.root"),
     "HNL_electronType_mHNL4p0_pl1000_rwctau2000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL4p0_pl1000/HeavyNeutralLepton_Tree.root"),
     ##"HNL_electronType_mHNL7p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL7p0_pl100/HeavyNeutralLepton_Tree.root"),
     #"HNL_electronType_mHNL7p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL7p0_pl1000/HeavyNeutralLepton_Tree.root"),
     ##"HNL_electronType_mHNL7p0_pl10000"  :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL7p0_pl10000/HeavyNeutralLepton_Tree.root"),
     ##"HNL_electronType_mHNL10p0_pl100"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL10p0_pl100/HeavyNeutralLepton_Tree.root"),
     #"HNL_electronType_mHNL10p0_pl1000"  :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL10p0_pl1000/HeavyNeutralLepton_Tree.root"),
     #"HNL_electronType_mHNL10p0_pl10000":glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL10p0_pl10000/HeavyNeutralLepton_Tree.root"),
     ##  
     'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8':glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/WJetsToLNu/HeavyNeutralLepton_Tree.root")       , 
     'EGamma_2018A': glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/EGamma_2018A/HeavyNeutralLepton_Tree_*.root"),
     'EGamma_2018B': glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/EGamma_2018B/HeavyNeutralLepton_Tree_*.root"),
     'EGamma_2018C': glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/EGamma_2018C/HeavyNeutralLepton_Tree_*.root"),
     'EGamma_2018D': glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/EGamma_2018D/HeavyNeutralLepton_Tree_*.root"),
}
signals_muon = {
      "HNL_muonType_mHNL1p0_pl10"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL1p0_pl10/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL1p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL1p0_pl100/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL1p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL1p0_pl1000/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL2p0_pl10"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL2p0_pl10/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL2p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL2p0_pl100/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL2p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL2p0_pl1000/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL4p0_pl10"     :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL4p0_pl10/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL4p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL4p0_pl100/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL4p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL4p0_pl1000/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL7p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL7p0_pl100/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL7p0_pl1000"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL7p0_pl1000/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL7p0_pl10000"  :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL7p0_pl10000/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL10p0_pl100"   :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL10p0_pl100/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL10p0_pl1000"  :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL10p0_pl1000/HeavyNeutralLepton_Tree_*.root"),
     "HNL_muonType_mHNL10p0_pl10000":glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_muonType_mHNL10p0_pl10000/HeavyNeutralLepton_Tree_*.root"),

}
data_muon = {
     'Muon_2018D': glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/SingleMuon_2018D/HeavyNeutralLepton_Tree_*.root"),

}
test = {
     "HNL_electronType_mHNL2p0_pl100"    :glob.glob("/eos/uscms/store/user/lpclonglived/HNL/skim/HNL_electronType_mHNL2p0_pl100/HeavyNeutralLepton_Tree_1.root"),
}


def writejson(fileset,fout):
    outf = open(os.path.expandvars(fout),"w")
    finaljson={}
    for sample, filelist in fileset.items():
        print("working on ",sample)
        newpaths = []
        for f in filelist:
            #newpaths.append(f.replace("/eos/uscms/",'root://xcache//')) ## coffea-casa
            newpaths.append(f.replace("/eos/uscms/",'root://cmseos.fnal.gov//'))
        finaljson[sample] = newpaths
    outf.write((json.dumps(finaljson,indent=4)))
    print("Json written to :", fout)
    return

# For coffea-casa
def replaceFilepath(inputjson):
    outputjson = inputjson.replace(".json","_casa.json")
    with open(inputjson) as f:
        d = json.load(f)
    output ={}
    for k,flist in d.items():
        output[k]=[f.replace("cmsxrootd.fnal.gov","xcache") for f in flist]
    with open(outputjson, 'w') as outfile:
        outfile.write(json.dumps(output,indent=4))

#writejson(test,"test.json")
#writejson(signals,"signals.json")
#writejson(signals,"signals_skim.json")
#writejson(signals_muon,"signals_muon.json")
writejson(data_muon,"data_muon.json")
#replacejson("signals.json")
