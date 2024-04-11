import glob
import sys, os 
from optparse import OptionParser,OptionGroup
import json
import numpy as np

def exec_me(command, dryRun=False):
    print(command)
    if not dryRun:
        os.system(command)

def write_condor(exe='runjob', files = [], dryRun=True):
    fname = '%s.jdl' % exe
    out = """universe = vanilla
executable = $(filename)
Should_Transfer_Files = YES 
WhenToTransferOutput = ON_EXIT_OR_EVICT
Transfer_Input_Files = $(filename),{files}
Output = {exe}.$(Process).$(Cluster).stdout
Error  = {exe}.$(Process).$(Cluster).stdout
Log    = {exe}.$(Process).$(Cluster).log
Queue filename matching {exe}_*.sh
""".format(exe=exe, files=','.join(files))
    
    with open(fname, 'w') as f:
        f.write(out)
    if not dryRun:
        os.system("condor_submit %s" % fname)


def write_bash(temp = 'runjob.sh',command = '', hadd="", files=[],eoscp=""):
    CMSSW="CMSSW_10_2_13"
    out = '#!/bin/bash\n'
    #out += 'while getopts e:m:t: flag\n'
    #out += 'do\n'
    #out += '    case "${flag}" in\n'
    #out += '        e) fe=${OPTARG};;\n'
    #out += '        m) fmu=${OPTARG};;\n'
    #out += '        t) ftau=${OPTARG};;\n'
    #out += '    esac\n'
    #out += 'done\n'
    out += 'date\n'
    out += 'MAINDIR=`pwd`\n'
    out += 'ls\n'
    out += 'voms-proxy-info --all\n'
    out += 'source /cvmfs/cms.cern.ch/cmsset_default.sh\n'
    out += 'export CWD=${PWD}\n'
    out += 'export PATH=${PATH}:/cvmfs/cms.cern.ch/common\n'
    out += 'export SCRAM_ARCH=slc7_amd64_gcc700\n'
    out += 'tar -xf CMSSW_10_2_13.tgz\n'
    for f in transfer_files:
        basename = os.path.basename(f)
        if ".py" in basename or ".txt" in basename:
            out += 'mv %s %s/src\n'%(basename,CMSSW) ## move executable to folder
    out += 'cd CMSSW_10_2_13/src\n'
    out += 'rm -r  CombineHarvester\n' ##remove combineHavester for now
    out += 'scramv1 b ProjectRename\n'
    out += 'eval `scramv1 runtime -sh` # cmsenv\n'
    out += 'echo $CMSSW_BASE "is the CMSSW we have on the local worker node"\n'
    out += 'echo "After Untar: "\n'
    out += 'ls -lrth\n'
    out += command + '\n'
    out += 'echo "Command: %s"\n'%command
    out += 'echo "hadd: %s"\n'%hadd
    out += hadd + '\n'
    out += 'cd ${CWD}\n'
    out += 'mv ./CMSSW_10_2_13/src/*.root .\n'        #collect output
    if eoscp!="":
        out += 'echo "coping to eos: "+%s  \n'%eoscp
        out +=  eoscp + '\n'
    out += 'echo "Inside $MAINDIR:"\n'
    out += 'ls\n'
    out += 'echo "DELETING..."\n'
    out += 'rm -rf %s\n'%CMSSW
    out += 'rm -rf *.pdf *.C core*\n'
    if eoscp!="":
        out += 'cd $MAINDIR  \n'
        out += 'echo "remove output local file"  \n'
        out += 'rm -rf *.root \n'

    out += 'ls\n'
    out += 'date\n'
    with open(temp, 'w') as f:
        f.write(out)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--clean', dest='clean', action='store_true',default = False, help='clean submission files', metavar='clean')
    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='write submission files only', metavar='dryRUn')
    parser.add_option('-i', '--idir', dest='idir', default='./', help='directory of input cards', metavar='idir')

    script_group  = OptionGroup(parser, "script options")
    parser.add_option_group(script_group)

    (options, args) = parser.parse_args()
    dryRun= options.dryRun 
    cardpath= options.idir
    #cardpath = "./combine/HNL_datacards/ele_v19/"

    outpath= cardpath+"hybrid/"

    allcards = glob.glob(cardpath+"*_comb.txt")
    ncards   = len(allcards)
    basedir  = os.path.dirname  
    print(" Total cards = ",ncards)

    if not os.path.exists(outpath):
      exec_me("mkdir -p %s"%(outpath), False)
    ##Small files used by the exe
    transfer_files = [
        os.getcwd()+"/hybrid.py",
        os.getcwd()+"/CMSSW_10_2_13.tgz"
    ]

    # transfer all cards for all jobs
    transfer_files += [ os.getcwd()+"/"+card for card in allcards ]
    # transfer all cards for runjob_%i.sh 
    #transfer_files += [ os.getcwd()+"/runjob_%s.sh"%i for i,card in enumerate(allcards)]

    os.chdir(outpath)
    print("submitting jobs from : ",os.getcwd())
    for i,card in enumerate(allcards):
        exe = "runjob"
        fname   = os.path.basename(card)
        SignalName    = fname.replace("_comb.txt","")
        command = f"python hybrid.py --card {fname}"
        hadd    =  f"hadd higgsCombine_{SignalName}_hybridnew_merged.root *.root"
        write_bash(exe+"_%i.sh"%i,  command, hadd, transfer_files,"")
    write_condor( exe,  transfer_files, dryRun)
