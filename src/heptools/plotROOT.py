import mplhep as hep
import matplotlib.pyplot as plt 
from ROOT import TH1

data_label = {"text"  : "Preliminary",
                "data": False,
                "lumi": 140,
                "year": 2016}

def atlas_1d(input_histograms, 
             norm:bool      = False,
             input_label    = data_label, 
             legend         = [], 
             **kwargs):

    """
    Script for generating (non-ratio) plot of ROOT TH1 histograms
    """

    # Style
    hep.style.use(hep.style.ATLAS) 

    if isinstance(input_histograms,TH1):
        list_of_histograms = [input_histograms]
    elif isinstance(input_histograms,list):
        list_of_histograms = input_histograms
    else:
        raise TypeError("The input_histograms field can only be a TH1 histogram or a list of TH1 histograms")

    if norm:
        hists_to_plot = []
        for h in list_of_histograms:
            hc = h.Clone(h.GetName())
            hc.Scale(1/hc.Integral())
            hists_to_plot.append(hc)
    else:
        hists_to_plot = list_of_histograms

    # Plot
    hep.histplot(hists_to_plot,**kwargs)

    data_label.update(input_label)
    hep.atlas.label(data_label["text"], data=data_label["data"], lumi=data_label["lumi"], year=data_label["year"])
        
    # Legend
    if isinstance(legend,bool) and legend==True:
        legend = [h.GetName() for h in list_of_histograms]

    plt.gca().legend(legend)

    return plt

