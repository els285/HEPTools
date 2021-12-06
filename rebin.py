'''
Ethan Simpson December 6th 2021
Rebinning two-dimensional ROOT histograms
Function requires: - original ROOT histogram - new axes 
Note that the new axes must a subset of the original axes.
'''

from ROOT import TH2D
import numpy as np

from functools import wraps
from time import process_time

# This turns a ROOT TAxis object into a numpy array
get_axis = lambda hist_axis: np.asarray([hist_axis.GetBinUpEdge(binn) for binn in range(0,hist_axis.GetNbins()+1)])

# Decorator for measuring execution time of the rebin2D function - not necessary
def measure_time(func):
    @wraps(func)
    def _time_it(*args, **kwargs):
        start = int(round(process_time() * 1000))
        try:
            return func(*args, **kwargs)
        finally:
            end_ = int(round(process_time() * 1000)) - start
            print(
                f"Total execution time {func.__name__}: {end_ if end_ > 0 else 0} ms"
            )

    return _time_it

# Main function
@measure_time
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

