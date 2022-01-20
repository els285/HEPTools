import julia
from julia import Main


import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

from heptools.eft.ParsefromJSON import parse_from_json
from heptools.eft.Plot_1D_Posterior import _1D_LHC_Posterior_Plot as LHC_p1D



def fetch_1d_posterior(filename):

    get_sample = Main.include("parse_sample_function2.jl")
    get_weights = Main.include("parse_weights.jl")

    # Requires the Wilson coefficients to be in the order they are in runEFT...jl

    list_of_wc = ["cHQ1","cHQ3","cHt"]

    data = {}
    hists = {}

    for i,wc in enumerate(list_of_wc):
        arr=np.asarray(get_sample(filename,i+1))
        data[wc]=arr
        hists[wc]=np.histogram(arr,bins=50)

    df = pd.DataFrame(data)

    df["weights"] = get_weights(filename)

    return df,hists


def make_1D_posterior_LHC_plot(WC,experiment,samples_filename,limits_filename):

    samples_df , hists = fetch_1d_posterior(samples_filename)

    limits_df = parse_from_json(limits_filename).T

    Plot_Object = LHC_p1D(samples_df[WC] , samples_df["weights"] , limits=limits_df[WC],experiment="ATLAS")
    Plot_Object.make_plot(include_global_mode=True,include_prior=True)
    Plot_Object.additional_label(WC+" operator")

    return Plot_Object



