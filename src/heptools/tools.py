import numpy as np

##############################################################
###---                 ROOT-based tools                 ---###
##############################################################

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


#################### Rebinning tools ######################

def rebin2D(hist2D,Xaxis,Yaxis):

    new_hist2D = TH2D(hist2D.GetName()+'_rebinned','',len(Xaxis)-1,Xaxis,len(Yaxis)-1,Yaxis)

    oXaxis = get_axis(hist2D.GetXaxis())
    oYaxis = get_axis(hist2D.GetYaxis())

    assert all([x in oXaxis for x in Xaxis]), "Values in the new X-axis are not contained in the original X-axis"
    assert all([y in oYaxis for y in Yaxis]), "Values in the new Y-axis are not contained in the original Y-axis"

    oXaxis_split = np.asarray([np.where(oXaxis==x)[0][0] for x in Xaxis])
    oYaxis_split = np.asarray([np.where(oYaxis==y)[0][0] for y in Yaxis])

    for i in range(0,len(oXaxis_split)-1):
        for j in range(0,len(oYaxis_split)-1):

            new_bin_total = 0
            for s in range(oXaxis_split[i],oXaxis_split[i+1]):
                for t in range(oYaxis_split[j],oYaxis_split[j+1]):
                    new_bin_total += hist2D.GetBinContent(s+1,t+1)
            
            new_hist2D.SetBinContent(i+1,j+1,new_bin_total)

    return new_hist2D


def rebin2D_Nbins(hist2D,XNbins,YNbins):

    oXaxis = get_axis(hist2D.GetXaxis())
    oYaxis = get_axis(hist2D.GetYaxis())

    Xaxis = np.linspace(oXaxis[0],oXaxis[-1],XNbins+1)
    Yaxis = np.linspace(oXaxis[0],oXaxis[-1],YNbins+1)

    return rebin2D(hist2D,Xaxis,Yaxis)

def rebin1D(hist1D,Xaxis):

    new_hist1D = TH1D(hist1D.GetName()+'_rebinned','', len(Xaxis)-1, Xaxis)

    oXaxis = get_axis(hist1D.GetXaxis())

    assert all([x in oXaxis for x in Xaxis]), "Values in the new X-axis are not contained in the original X-axis"

    oXaxis_split = np.asarray([np.where(oXaxis==x)[0][0] for x in Xaxis])

    for i in range(0,len(oXaxis_split)-1):
        
        new_bin_total = 0
        for s in range(oXaxis_split[i],oXaxis_split[i+1]):
            new_bin_total += hist1D.GetBinContent(s+1)
        
        new_hist1D.SetBinContent(i+1,new_bin_total)

    return new_hist1D


def rebin1D_Nbins(hist1D,Nbins):

    oXaxis = get_axis(hist1D.GetXaxis())

    Xaxis = np.linspace(oXaxis[0],oXaxis[-1],Nbins+1)

    return rebin1D(hist1D,Xaxis)




##############################################################
###---               scikit-HEP-based tools             ---###
##############################################################

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
    column_names = {df_out.columns[og]:new for og,new in enumerate([t[1] for t in inputs])}
    df_final = df_out.rename(columns=column_names)#, inplace=True)
    return df_out


