from optparse import OptionParser
import glob,os

def Run(cmd,test=False):
    print(cmd)
    if not test: os.system(cmd)

def comb(card,test):

    basedir = os.path.dirname(card)
    name = card.split("/")[-1].replace(".txt","")

    outdir = basedir
    quantiles       = [-1,0.16,0.5,0.84,0.975,0.025]

    for q in quantiles:
        cmd = "combine -M HybridNew --LHCmode LHC-limits {card} -n _{name}_hybridnew --setParameters norm=1 --freezeParameter norm --saveHybridResult".format(name=name,card=card)
        if not q==-1:
            cmd += " --expectedFromGrid {q} --fork 16".format(q=q)
        Run(cmd,test)
    #Run("mv higgsCombine_%s_hybridnew.HybridNew.mH120*.root %s"%(name,outdir),test)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='dryRun')
    parser.add_option('--card', dest='card', default = "card.txt", help='input combine card')

    (options, args) = parser.parse_args()

    #comb("./combine/HNL_datacards/muon_v25/HNL_muonType_mHNL2p0_pl1000_comb.txt",options.dryRun)
    #comb("./combine/HNL_datacards/ele_v19/HNL_electronType_mHNL2p0_pl1000_comb.txt",options.dryRun)
    comb(options.card,options.dryRun)
