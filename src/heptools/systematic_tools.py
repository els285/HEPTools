# Systematics Plotting tools

from ROOT import TCanvas, TFile 
import os 

from os import listdir
from os.path import isfile, join


import matplotlib.pyplot as plt

from heptools.histplot.PyHist_Class import Histogram_Wrapper as Hist
from heptools.histplot.Plotting import Ratio_Plot_ROOT



def generate_symmetrised_hist(hist,NOMINAL_hist):

    """
    ROOT-based symmetrisation of a particular histogram around a Nominal histogram
    """

    SYM_hist = hist.Clone(hist.GetName()+"_symmetrised")

    for i in range(0,hist.GetNbinsX()):
        diff = hist.GetBinContent(i+1) - NOMINAL_hist.GetBinContent(i+1)
        SYM_hist.SetBinContent(i+1,NOMINAL_hist.GetBinContent(i+1)-diff)

    return SYM_hist




def plot_up_down_variation(NOMINAL_hist,UP_hist,DOWN_hist,**kwargs):


    """
    Generate an ATLAS "Up-Down" variation plot
    """

    x_nom = Hist(NOMINAL_hist,name="nominal",colour="black"         ,legend_entry="Nominal")
    x_up  = Hist(UP_hist,name="UP variation",colour="red"           ,legend_entry="UP")
    x_down = Hist(DOWN_hist,name="DOWN variation",colour="blue"     ,legend_entry="DOWN")

    norm=kwargs["norm"] if "norm" in kwargs else False

    p = Ratio_Plot_ROOT("A Plot",list_of_histograms=[x_nom,x_up,x_down],divisor=x_nom,normalise=norm)
    p.Initialise_Plot_Design("ATLAS")

    # plt.figure(figsize=(10,10)) 
    plt,ax,rax = p.Make_Ratio_Plot("basic-line")
    plt.gcf().set_size_inches(10, 10)
    p.Add_ATLAS_Label("Internal")
    # rax.set_ylim([0.99,1.01])
    ax.legend(loc=4)
    rax.set_xlabel(r"$\cos \phi$")

    return plt,ax,rax

