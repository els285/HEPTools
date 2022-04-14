import matplotlib.pyplot as plt

from heptools.histplot import PyHist_Class_2D
from heptools.histplot.PyHist_Class import Histogram_Wrapper
from heptools.histplot import RatioPlot 
from heptools.histplot import SinglePlot


def quick2D(ROOT_histogram,**kwargs):

    pyhist = PyHist_Class_2D.Histogram_Wrapper(ROOT_histogram,"Name")
    fig,ax,pc=pyhist.plot_2d(normed=True,**kwargs)

    return pyhist,fig,ax,pc


def ATLAS_1D(list_of_ROOT_histograms,**kwargs):

    """ Pass a list of ROOT histograms and generate a plot"""
    import matplotlib.pyplot as plt
    default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']    

    list_of_histograms = []
    for ROOThist,col in zip(list_of_ROOT_histograms,default_colours):
        x=Histogram_Wrapper(ROOThist, name = ROOThist.GetName()  ,colour=col , legend_entry = ROOThist.GetName()) 
        list_of_histograms.append(x)

    p,plt = SinglePlot.standard_ATLAS_plot(list_of_histograms,**kwargs)

    return p,plt


def ATLAS_ratio_1D(list_of_ROOT_histograms,**kwargs):

    """ Pass a list of ROOT histograms and generate a ratio plot"""

    import matplotlib.pyplot as plt
    default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']    

    list_of_histograms = []
    for ROOThist,col in zip(list_of_ROOT_histograms,default_colours):
        x=Histogram_Wrapper(ROOThist, name = ROOThist.GetName()  ,colour=col , legend_entry = ROOThist.GetName()) 
        list_of_histograms.append(x)

    p,plt = RatioPlot.standard_ATLAS_ratio_plot(list_of_histograms,**kwargs)

    return p,plt
