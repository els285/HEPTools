import mplhep as hep
import matplotlib.pyplot as plt 

def atlas_1d(list_of_hists,**kwargs):

    # Main
    hep.style.use(hep.style.ATLAS) 
    hep.histplot(list_of_hists)

    # Accompanying text
    if "text" in kwargs:
        hep.atlas.text(kwargs["text"])
    else:
        hep.atlas.text("Preliminary")


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

