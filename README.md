# HEPTools

A collection of tools for high-energy particle physics analyses.

**Crucial: ROOT is required to be installed and compatible with Python3. I only have ROOT built for Python3.6**


## Installation

Install via
```
python3.6 -m pip install git+https://github.com/ethansimpson285/HEPTools
```


What follows is a descrption of the various tools...


## Plotting Tools

* Functionality to plot a one-dimensional histograms in a variety of designs, including LHC-experiment-based as facilitated through the `mplhep` module`
* Functionality to plot two-dimensional histograms
* Functionality for generating various effective field theory interpretation plots.
(All of the above are matplotlib based, with interest in developing online tools - likely through `plotly` - down-the-line.

### Single Histogram Plotting
This is designed to be as simple as possible. Do the following to plot a ROOT histogram:
```python3
from heptools.quick_plot import atlas_1D_plot 
p,plt = atlas_1D_plot(ROOT_hist)
```
Pass a list of histograms to plot them on top of one another: `atlas_1D_plot([hist1,hist2,hist3])`.

For greater control, you can generate the `Histogram_Wrapper` objects yourself before plotting.

### Ratio Plots
The analogous, simple ratio plotting script is:
```python3
from heptools.quick_plot import atlas_1D_ratio_plot
p,plt = atlas_1D_ratio_plot([hist1,hist2]),divisor=hist1)
```
If no `divisor` is specified, the first histogram in the list is taken as the divisor.


## Parsers

* `LHE2ROOT` is a parser designed to convert LHE (Les Houches accord) file types, common outputs of HEP Monte Carlo event generators, into the ROOT file format.

## Skimming
For trimming branches of TTrees.

* `generic_tree_skim` applies a basic skimming based on a list of branches passed as an argument:
  ```python3
    _file0 = TFile("~/mc15_13TeV.410011.PwPyEG_P2012_singletop_tchan_lept_top.1lep_raw.root")
    tree = _file0.Get("mini")
    branches2keep = ["eventNumber","lep_type","lep_pt"]
    generic_tree_skim(tree,branches2keep,cut_off=200000,file_name="hi.root",tree_name="parton_tree")
  ```
