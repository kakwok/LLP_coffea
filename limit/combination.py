import os
import sys
from optparse import OptionParser


#input file
def add_prefix(cardfile,nuisances,prefix):
    fin = open(cardfile, "rt")
    fout = open(cardfile.replace(".txt","renamed.txt"), "wt")
    for line in fin:
    	temp = line
    	for nuisance in nuisances:
    		temp=temp.replace(nuisance, prefix+nuisance)
    	fout.write(temp)
    fin.close()
    fout.close()
    return

if __name__ == "__main__":

    parser = OptionParser(
        usage="python combination prefix1=[card1] prefix2=[csc_card] [output_card]",
    )
    (options, args) = parser.parse_args()

    print(args)
    ##list of nuisances to be renamed
    dt_nuisances = ['clusterEff_unc','cut_based_eff_unc','jetVeto','muonVeto','rechitVeto','time','time_spread','readout','NA','NB','NC','ND']
    csc_nuisances = ['clusterEff_unc','cut_based_eff_unc','jetVeto','muonVeto','rechitVeto','time','time_spread','readout','NA','NB','NC','ND']

    # correlated: higgs xsec, pdf, shape, lumi, JEC, MC stats
    # uncorrelated: CSC hit related, add csc/dt prefix
    
    ##
    ## correlated   : jetVeto (to be checked)
    ## uncorrelated : muonVeto? time, time_spread
    
    ## Dark shower signal
    ## Double counting: overlap region and 2 clusters in DT/CSC

    if not args:
        raise(RuntimeError, "No input datacards specified.")

    (prefix1,dt_card)  = args[0].split("=")
    (prefix2,csc_card) = args[1].split("=")
    output_card        = args[2]

    add_prefix(dt_card,dt_nuisances,prefix1+"_")
    add_prefix(csc_card,csc_nuisances,prefix2+"_")

    print("Writing output datacard:",output_card)
    os.system('combineCards.py {p1}={dc1} {p2}={dc2}  > temp_output.txt'.format(
                p1=prefix1,p2=prefix2,
                dc1=dt_card.replace(".txt","renamed.txt"),
                dc2=csc_card.replace(".txt","renamed.txt")
                ))
    with open('temp_output.txt', 'r') as read_obj, open(output_card, 'w') as write_obj:
        write_obj.write(open(csc_card, "r").readline() + '\n')
        for line in read_obj:
            write_obj.write(line)
    
    os.system('rm {}'.format(dt_card.replace(".txt","renamed.txt")))
    os.system('rm {}'.format(csc_card.replace(".txt","renamed.txt")))
    os.system('rm temp_output.txt')
