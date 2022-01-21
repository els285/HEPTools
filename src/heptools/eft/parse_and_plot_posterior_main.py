# Skip the julia code altogether
import julia
from julia import Main

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

from Fetch1DPosterior import fetch_1d_posterior

import time 

st = time.time()





def pull_data(WC,samples_filename,limits_filename):

    posterior_df , hists = fetch_1d_posterior(samples_filename)

    limits_df = 

    
        

df,hists = fetch_1d_posterior("../results_ptz_laurynas_ethantest/BATobjects/samples.jld2")

print("df made")
print(time.time()-st)

from Plot_1D_Posterior import _1D_LHC_Posterior_Plot as LHC_p1D

plot = LHC_p1D(df.cHt,limits={"1sigma":(-2.4,6.4),"2sigma":(-5.8,9.2)},experiment="ATLAS")
plot.make_multi_plot()
plot.include_metadata(text="",data=True,lumi=50,year=2017)
plot.additional_label("cHt operator")


plot.ax.set_xlim(-50,50)
plt.show()
input()
