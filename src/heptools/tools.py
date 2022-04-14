import numpy as np

def bin_values(hist):
    return np.asarray([hist.GetBinContent(binn+1) for binn in range(0,hist.GetNbinsX())])

def bin_errors(hist):
    return np.asarray([hist.GetBinError(binn+1) for binn in range(0,hist.GetNbinsX())])

def bin_edges(hist):
    return np.asarray([hist.GetXaxis().GetBinUpEdge(binn) for binn in range(0,hist.GetXaxis().GetNbins()+1)])

def bin_centres(hist):
    return np.asarray([hist.GetXaxis().GetBinCenter(binn+1) for binn in range(0,hist.GetNbinsX())])


def compute_ratio(ROOT_hist_numer,ROOT_hist_denom):

    d_hist = ROOT_hist_numer.Clone()
    d_hist.Divide(ROOT_hist_numer,ROOT_hist_denom)
    return d_hist


def key_names(file):
    return [k.GetName() for k in file.GetListOfKeys()]

def branch_names(ttree):
    return [k.GetName() for k in ttree.GetListOfBranches()]        