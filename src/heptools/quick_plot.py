import matplotlib.pyplot as plt
from ROOT import TH1

from heptools.histplot import PyHist_Class_2D
from heptools.histplot.PyHist_Class import Histogram_Wrapper
from heptools.histplot import RatioPlot 
from heptools.histplot import SinglePlot
from heptools.histplot import DataMCPlot




def quick2D(ROOT_histogram,**kwargs):

    pyhist = PyHist_Class_2D.Histogram_Wrapper(ROOT_histogram,"Name")
    fig,ax,pc=pyhist.plot_2d(normed=True,**kwargs)

    return pyhist,fig,ax,pc


def atlas_1D_plot(input_histograms,**kwargs):

    """ Pass a list of ROOT histograms and generate a plot"""

    import matplotlib.pyplot as plt
    default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']    

    if isinstance(input_histograms,TH1):
        list_of_ROOT_histograms = [input_histograms]
    elif isinstance(input_histograms,list):
        list_of_ROOT_histograms = input_histograms
    else:
        print("Histogram input must be Histogram_Wrapper class or list of Histogram_Wrappers")

    list_of_histograms = []
    for ROOThist,col in zip(list_of_ROOT_histograms,default_colours):
        x=Histogram_Wrapper(ROOThist, name = ROOThist.GetName()  ,colour=col , legend_entry = ROOThist.GetName()) 
        list_of_histograms.append(x)

    p,plt = SinglePlot.standard_ATLAS_plot(list_of_histograms,**kwargs)

    return p,plt


def atlas_1D_ratio_plot(list_of_ROOT_histograms,**kwargs):

    """ Pass a list of ROOT histograms and generate a ratio plot"""

    import matplotlib.pyplot as plt
    default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']    

    list_of_histograms = []
    for ROOThist,col in zip(list_of_ROOT_histograms,default_colours):
        x=Histogram_Wrapper(ROOThist, name = ROOThist.GetName()  ,colour=col , legend_entry = ROOThist.GetName()) 
        list_of_histograms.append(x)

    p,plt = RatioPlot.standard_ATLAS_ratio_plot(list_of_histograms,**kwargs)

    return p,plt



def atlas_dataMC_plot(list_of_MC_ROOT_histograms,data_ROOT_histogram):

    import matplotlib.pyplot as plt
    default_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']    

    list_of_MC_histograms = []
    for ROOThist,col in zip(list_of_MC_ROOT_histograms,default_colours):
        x=Histogram_Wrapper(ROOThist, name = ROOThist.GetName()  ,colour=col , legend_entry = ROOThist.GetName()) 
        list_of_MC_histograms.append(x)   

    dataPyHist = Histogram_Wrapper(data_ROOT_histogram)

    p,plt = DataMCPlot.standard_ATLAS_dataMC_plot(list_of_MC_histograms,dataPyHist)

    plt.show()
    input()


