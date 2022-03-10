# Running expected limit

To make the datacards, we need the following inputs:
 1. Selection cuts on ABCD plane
    - Selection cuts are done with this notebook:`notebook/ClosureChecksAndCutOptimization.ipynb`
 2. Background prediction (A,B,C,D)
    - Background prediction are done once per channel with the optimized cuts using the function [predIntimeFromOOT](https://github.com/kakwok/LLP_coffea/blob/refactor/limit/writeABCD.py#L201)
    - Background values for muon channels are hard-coded [here](https://github.com/kakwok/LLP_coffea/blob/refactor/limit/writeABCD.py#L292-L293)    
 3. Signal prediction (A,B,C,D)
    - Repeat the process for all signal samples using the same cuts

## Signal predictions

After preparing the signal histograms with the processor, run a **coffea env**
```
python writeABCD.py --loadhist
```
This prints out the yields for each signal samples to the terminal, which can be copied to the [calSigal](https://github.com/kakwok/LLP_coffea/blob/refactor/limit/writeABCD.py#L66) function.

Run the script in a separate terminal **without coffea env**, to write all the data cards
```
python writeABCD.py
```

The manual copy is done for now because combine is not compatible with coffea environment.
We could improve remove the manual copy by writing the the signal yields to a JSON file and read by the script as input next.  

## Limit plotting

This is done by the [limitPlotting ](https://github.com/kakwok/LLP_coffea/blob/refactor/notebook/LimitPloting.ipynb) notebook

The key functions are the following:
 - `loadLimitFiles`: Read combine card limits into numpy arrays
    - Called once per-mass point
 - `plotlimit` : Make brazil plot given the combine-style numpy arrays
 - `interpolate1D`: Returns limit points, ctau, and the interpolated function to evaluate the limit in other ctau points
 - `intersect`: Find the intersection points between theory/ctau point
 -  Helper functions:
    - `f_xsec(m)`, `f_ctau(m)` : returns a function for xsec/ctau given HNL mass `m` 
  
