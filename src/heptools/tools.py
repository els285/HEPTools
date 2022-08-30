import numpy as np

def bin_values(hist):
    return np.asarray([hist.GetBinContent(binn+1) for binn in range(0,hist.GetNbinsX())])

def bin_errors(hist):
    return np.asarray([hist.GetBinError(binn+1) for binn in range(0,hist.GetNbinsX())])

def bin_edges(hist):
    return np.asarray([hist.GetXaxis().GetBinUpEdge(binn) for binn in range(0,hist.GetXaxis().GetNbins()+1)])

def bin_centres(hist):
    return np.asarray([hist.GetXaxis().GetBinCenter(binn+1) for binn in range(0,hist.GetNbinsX())])

def normalise(hist):
    """ Takes a ROOT histogram and normalised it"""
    h1dN = hist.Clone(hist.GetName()+"_norm")
    h1dN.Scale(1/hist.Integral())
    return h1dN    


def compute_ratio(ROOT_hist_numer,ROOT_hist_denom):

    d_hist = ROOT_hist_numer.Clone()
    d_hist.Divide(ROOT_hist_numer,ROOT_hist_denom)
    return d_hist


def key_names(file):
    return [k.GetName() for k in file.GetListOfKeys()]

def branch_names(ttree):
    return [k.GetName() for k in ttree.GetListOfBranches()]        


def merge_dataframes(dfs,common_key: str):
    import functools 
    import pandas as pd
    return functools.reduce(lambda left, right: pd.merge(left, right, on=common_key), dfs)


def friend_ROOT_trees(inputs,**kwargs):

    """
    Takes input of form: [(file1,tree1),(file2,tree2),...]
    Collects all matched events in dataframe - default collection key is "eventNumber"
    """

    import uproot
    import awkward as ak
    import functools
    import pandas as pd

    common_key = kwargs["common_key"] if "common_key" in kwargs else "eventNumber"

    dfs = []

    for file in inputs:
        filename = file[0]
        treename = file[1]

        # Extract the eventNumber TBranch as a pandas DataFrame
        with uproot.open(f"{filename}:{treename}") as events:
            eventNumbers = events.arrays([common_key],library="pd")
        eN = eventNumbers.reset_index(level=1, drop=True)

        # Extract the full desired Branches as awkward array, and parse
        with uproot.open(f"{filename}:{treename}") as events:
            list_of_branches = kwargs["list_of_branches"] if "list_of_branches" in kwargs else events.keys()
            outputs = events.arrays(list_of_branches,library="ak")
        arrow_array = ak.to_arrow(outputs)
        df = arrow_array.to_pandas()
        dff = df.to_frame()

        # Merge and set eventNumber as index
        merged = pd.concat([eN,dff],axis=1)
        final = merged.set_index([common_key])

        dfs.append(final)

    df_out = functools.reduce(lambda left, right: pd.merge(left, right, how="inner",on=common_key), dfs)
    return df_out


