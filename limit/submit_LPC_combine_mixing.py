import glob
import sys, os 
from optparse import OptionParser,OptionGroup
import json
import numpy as np


def exec_me(command, dryRun=False):
    print(command)
    if not dryRun:
        os.system(command)

def write_condor(njobs, exe='runjob', files = [],mixing="0.1 0.5 0.4", dryRun=True):
    fname = '%s.jdl' % exe
    out = """universe = vanilla
Executable = {exe}.sh
Should_Transfer_Files = YES 
WhenToTransferOutput = ON_EXIT_OR_EVICT
Transfer_Input_Files = {exe}.sh,{files}
Output = {exe}.$(Process).$(Cluster).stdout
Error  = {exe}.$(Process).$(Cluster).stdout
Log    = {exe}.$(Process).$(Cluster).log
Arguments = -e $(fe) -m $(fmu) -t $(ftau)
Queue fe,fmu,ftau from(
  {mixing}
)""".format(exe=exe, files=','.join(files), njobs=njobs, mixing=mixing)
    
    with open(fname, 'w') as f:
        f.write(out)
    if not dryRun:
        os.system("condor_submit %s" % fname)


def write_bash(temp = 'runjob.sh',  command = '', files=[],eoscp=""):
    CMSSW="CMSSW_10_2_13"
    out = '#!/bin/bash\n'
    out += 'while getopts e:m:t: flag\n'
    out += 'do\n'
    out += '    case "${flag}" in\n'
    out += '        e) fe=${OPTARG};;\n'
    out += '        m) fmu=${OPTARG};;\n'
    out += '        t) ftau=${OPTARG};;\n'
    out += '    esac\n'
    out += 'done\n'
    out += 'date\n'
    out += 'MAINDIR=`pwd`\n'
    out += 'ls\n'
    out += 'echo ${flag}\n'
    out += 'echo ${fe}\n'
    out += 'voms-proxy-info --all\n'
    out += 'source /cvmfs/cms.cern.ch/cmsset_default.sh\n'
    out += 'export CWD=${PWD}\n'
    out += 'export PATH=${PATH}:/cvmfs/cms.cern.ch/common\n'
    out += 'export SCRAM_ARCH=slc7_amd64_gcc700\n'
    out += 'tar -xf CMSSW_10_2_13.tgz\n'
    for f in transfer_files:
        basename = os.path.basename(f)
        if ".py" in basename or ".json" in basename:
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
    out += 'cd ${CWD}\n'
    out += 'mv ./CMSSW_10_2_13/src/HNL_datacards_MuE/*.root .\n'        #collect output
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
    parser.add_option('-o', '--odir', dest='odir', default='./', help='directory to write histograms/job output', metavar='odir')
    parser.add_option( '--njobs', dest='njobs', default=1, type="int", help='Number of jobs to split into', metavar='njobs')
    parser.add_option( '--exe', dest='exe', default="Runllp_hnl_analyzer",  help='Executable name to run', metavar='exe')

    script_group  = OptionGroup(parser, "script options")
    script_group.add_option('--eleJson', dest='eleJson', default = "", help='json files for ele channels')
    script_group.add_option('--muJson', dest='muJson', default = "", help='json files for muon channels')
    parser.add_option('--testing', dest='testing', action='store_true',default = False, help='test single mixing only', metavar='testing')
    parser.add_option_group(script_group)

    (options, args) = parser.parse_args()
    dryRun= options.dryRun 

    maxJobs = options.njobs 

    outpath= options.odir

    ##Small files used by the exe
    transfer_files = [
        os.getcwd()+"/"+options.eleJson,
        os.getcwd()+"/"+options.muJson,
        os.getcwd()+"/writeABCD.py",
        os.getcwd()+"/combination.py",
        os.getcwd()+"/writeABCD_mixing_condor.py",
        os.getcwd()+"/CMSSW_10_2_13.tgz"
    ]

    if not os.path.exists(outpath):
      exec_me("mkdir -p %s"%(outpath), False)
    os.chdir(outpath)
    print("submitting jobs from : ",os.getcwd())

    eleFileName = os.path.basename(options.eleJson)
    muFileName = os.path.basename(options.muJson)
    ## JS TODO: Loop over phase space
    if options.testing:
        command = "python writeABCD_mixing_condor.py --eleJson {} --muJson {}  --fe $fe --fmu $fmu --ftau $ftau".format(eleFileName,muFileName)
         
        print("command to run: ", command)

        exe = "runjob"
        write_bash(exe+".sh",  command, transfer_files,"")
        write_condor(maxJobs, exe,  transfer_files,"0.5 0.4 0.1\n 0.4 0.2 0.4", dryRun)
    else:
        string_mixing_points = ""
        mixing_points = np.arange(0,1.02,0.02)
        for fe in mixing_points:

            for fmu in mixing_points:



                for ftau in mixing_points:

                    sum_of_floats = fe+fmu+ftau
                    sum_of_floats = round(sum_of_floats,2)

                    if sum_of_floats!=1: continue
                    string_mixing_points += str(round(fe,2)) +" " + str(round(fmu,2)) +" " +str(round(ftau,2)) +"\n"             
                    
                    
        command = "python writeABCD_mixing_condor.py --eleJson {} --muJson {} --fe $fe --fmu $fmu --ftau $ftau".format(eleFileName,muFileName)

        print("command to run: ", command)

        exe = "runjob"
        write_bash(exe+".sh",  command, transfer_files,"")
        write_condor(maxJobs, exe,  transfer_files,string_mixing_points, dryRun)

