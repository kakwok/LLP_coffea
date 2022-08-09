from coffea import lumi_tools

def printLumiByEras(lumivalues,eras):
    runs = lumivalues._lumidata[:, 0].astype('u4')
    lumis = lumivalues._lumidata[:, 1].astype('u4')
    for era in eras:
        runMask = (runs>=era["FromRun"]) & (runs<=era["ToRun"])
        if len(runs[runMask])>0:
            ll  = lumi_tools.LumiList(runs[runMask], lumis[runMask])
            print(era["Era"], "%.3f /fb"%(lumivalues.get_lumi(ll)/1000.))
        else:
            print(" No runs in %s from run %s to run %s"%(era["Era"],era["FromRun"],era["ToRun"]))
    print("Total  = %.3f /fb"%(lumivalues.get_lumi(lumi_tools.LumiList(runs,lumis))/1000.))
    return 


if __name__ == '__main__':

    lumivalues2016 = lumi_tools.LumiData("metadata/lumi2016_new.csv")
    lumivalues2017 = lumi_tools.LumiData("metadata/lumi2017.csv.gz")
    lumivalues2018 = lumi_tools.LumiData("metadata/lumi2018.csv")
    
    #from : https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis#DATA
    
    Eras2018 = [
        {"Era":"2018A","FromRun":315252,"ToRun":316995},
        {"Era":"2018B","FromRun":317080,"ToRun":319310},
        {"Era":"2018C","FromRun":319337,"ToRun":320065},
        {"Era":"2018D","FromRun":320673,"ToRun":325175},
     ]
    Eras2017 = [
        {"Era":"2017A","FromRun":294927,"ToRun":297019},
        {"Era":"2017B","FromRun":297046,"ToRun":299329},
        {"Era":"2017C","FromRun":299368,"ToRun":302029},
        {"Era":"2017D","FromRun":302030,"ToRun":303434},
        {"Era":"2017E","FromRun":303824,"ToRun":304797},
        {"Era":"2017F","FromRun":305040,"ToRun":306462},
     ]
    Eras2016 = [
        {"Era":"2016A","FromRun":271036,"ToRun":271658},
        {"Era":"2016B","FromRun":272007,"ToRun":275376},
        {"Era":"2016C","FromRun":275657,"ToRun":276283},
        {"Era":"2016D","FromRun":276315,"ToRun":276811},
        {"Era":"2016E","FromRun":276831,"ToRun":277420},
        {"Era":"2016F","FromRun":277772,"ToRun":278808},
        {"Era":"2016G","FromRun":278820,"ToRun":280385},
        {"Era":"2016H","FromRun":280919,"ToRun":284044},
     ]

    printLumiByEras(lumivalues2018,Eras2018)
    printLumiByEras(lumivalues2017,Eras2017)
    printLumiByEras(lumivalues2016,Eras2016)
