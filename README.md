# LLP_coffea
coffea analysis of LLP

## Getting coffea

For lpc condor, follow instructions at [lpcjobqueue](https://github.com/CoffeaTeam/lpcjobqueue)

```
voms-proxy-init
./shell                          ## enter virtual environment
[python makezip.py]              ## if you have updated anything in the `HNLprocessor` package
Singularity> python runHNL.py    ## switch to runLPC in the script 
```
Processing all samples takes about ~40 mins.

## Preparing filepaths

MC samples are located in the `store/group/lpclonglived/HNL` area, which is accessible via xrootd

To make the input file for the HNLprocessor, create the json file with
```
python filepaths.py
```

## Making corrections
```
python compile_corrections.py
```
