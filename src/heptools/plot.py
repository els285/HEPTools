import mplhep as hep
import matplotlib.pyplot as plt 

def atlas_1d(list_of_hists,**kwargs):

    """
    Script for generating (non-ratio) plot of ROOT TH1 histograms
    """

    # Style
    hep.style.use(hep.style.ATLAS) 

    if any([x in kwargs and x for x in ["norm","Norm","normalise","Normalise"]]):
        hists_to_plot = []
        for h in list_of_hists:
            hc = h.Clone(h.GetName())
            hc.Scale(1/hc.Integral())
            hists_to_plot.append(hc)
    else:
        hists_to_plot = list_of_hists

    # Plot
    hep.histplot(hists_to_plot)

    # Accompanying text
    if "text" in kwargs and not "label" in kwargs:
        hep.atlas.text(kwargs["text"])
    elif not "text" in kwargs and not "label" in kwargs:
        hep.atlas.text("Preliminary")

    # Label
    if "label" in kwargs:
        if isinstance(kwargs["label"],dict):
            data_label = kwargs["label"]
        elif isinstance(kwargs["label"],bool) and kwargs["label"]:
            data_label = {"text"  : "",
                            "data": False,
                            "lumi": 139,
                            "year": 2016}

        hep.atlas.label(data_label["text"], data=data_label["data"], lumi=data_label["lumi"], year=data_label["year"])

    # Legend
    if "legend" in kwargs:
        if isinstance(kwargs["legend"],list):
            legend_entries = kwargs["legend"]
        elif isinstance(kwargs["legend"],bool):
            if kwargs["legend"]:
                legend_entries = [h.GetName() for h in list_of_hists]
            else:
                legend_entries = []
    else:
        legend_entries = []

    plt.gca().legend(legend_entries)

    return plt

