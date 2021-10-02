# LLP_coffea
coffea analysis of LLP

## Getting coffea

### Interactive usage on LPC

First time setup:
```
tar -xvf coffeaenv.tar.gz
sed -i "3a export PYTHONPATH=${PWD}/coffeaenv/lib/python3.7/site-packages:\$PYTHONPATH" coffeaenv/bin/activate
source coffeaenv/bin/activate
```
After first time set-up, enable coffea env with 
```
source coffeaenv/bin/activate
``` 
#### Jupyter notebook access
You can launch Jupyter notebook server on LPC node and forward the port to your computer.

To do so, you need to connect to LPC with 
```
ssh -Y -L localhost:8887:localhost:8887 USERNAME@cmslpc-sl7.fnal.gov
``
After enabling coffea env, you can start the server with
```
jupyter notebook --no-browser --port=8887 --ip 127.0.0.1
```

## Preparing filepaths

MC samples are located in the `store/group/lpclonglived/HNL` area, which is accessible via xrootd

To make the input file for the HNLprocessor, create the json file with
```
python filepaths.py
```

## Making corrections

Correction files includes WpT reweighting, pileup corrections, and MC cross section.
```
python compile_corrections.py
```
## Running the processor

Edit `runHNL.py` with the output filename and the fileset json (made with `filepaths.py`)
For local testing with a small fileset:
```
python runHNL.py --test
```
For local testing with a fraction of the full fileset:
```
python runHNL.py --local
```

### Running with condor
To run with lpc condor, first follow instructions at [lpcjobqueue!!](https://github.com/CoffeaTeam/lpcjobqueue)

```
voms-proxy-init
./shell                          ## enter virtual environment
[python makezip.py]              ## if you have updated anything in the `HNLprocessor` package
Singularity> python runHNL.py    ## switch to runLPC in the script 
```
Processing all samples takes about ~10 mins.
