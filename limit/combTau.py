from optparse import OptionParser
import glob,os

def Run(cmd,test=False):
    print(cmd)
    if not test: os.system(cmd)

def comb(outdir,test):
    ele_odir = outdir+"ele/"
    mu_odir = outdir+"mu/"
   
    card_paths = glob.glob(ele_odir+"HNL*_comb.txt")
    for path in card_paths:
        name = path.split("/")[-1].replace("_comb.txt","")
        ##if not ("2p1" in name or "2p2" in name):  continue
        #if not ("Dirac" in name ):  continue
        #if not (("1p9" in name) or ("2p05" in name) or ("2p15" in name)): continue
        cmd = "python combination.py ele={ele_odir}{name}_comb.txt mu={mu_odir}{name}_comb.txt {odir}{name}_comb.txt".format(name=name,odir=outdir,ele_odir=ele_odir,mu_odir=mu_odir)
        Run(cmd,test)
        #cmd = "combine -M AsymptoticLimits {odir}{name}_comb.txt -n _{name}_comb --setParameters norm={norm} --freezeParameter norm -t -1 --toysFreq".format(name=name,odir=outdir,norm=1)
        cmd = "combine -M AsymptoticLimits {odir}{name}_comb.txt -n _{name}_comb --setParameters norm={norm} --freezeParameter norm ".format(name=name,odir=outdir,norm=1)
        Run(cmd,test)
        Run("mv higgsCombine_%s.AsymptoticLimits.mH120.root %s"%(name+"_comb",outdir),test)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='dryRun')

    (options, args) = parser.parse_args()

    #comb("./combine/HNL_datacards/tau_v1/",options.dryRun)
    #comb("./combine/HNL_datacards/tau_v2/",options.dryRun)
    #comb("./combine/HNL_datacards/tau_v3/",options.dryRun)
    #comb("./combine/HNL_datacards/tau_v4/",options.dryRun)
    #comb("./combine/HNL_datacards/tau_v5/",options.dryRun)
    #comb("./combine/HNL_datacards/tau_v6/",options.dryRun)
    #comb("./combine/HNL_datacards/tau_v7/",options.dryRun)
    comb("./combine/HNL_datacards/tau_v8/",options.dryRun)
