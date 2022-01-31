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


## Parsers

* `LHE2ROOT` is a parser designed to convert LHE (Les Houches accord) file types, common outputs of HEP Monte Carlo event generators, into the ROOT file format.

