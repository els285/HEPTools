
import matplotlib.pyplot as plt
from heptools.histplot.Plotting import HEP_Plot
from heptools.histplot.PyHist_Class import Histogram_Wrapper


class Single_Plot(HEP_Plot):
    def __init__(self,plot_title,**kwargs):
        super().__init__(plot_title,**kwargs) 


    def make_plot(self,plot_type):
        plotting_function = self.select_plot_type(plot_type)
        fig, ax = plt.subplots()
        self.axes = [ax]

        for HW in self.list_of_histograms:
            PH = HW.Norm_PyWrap_Hist if self.Normalised else HW.UnNorm_PyWrap_Hist
            plotting_function(ax,PH)
            self.do_legend()
        
        return plt, fig, ax


def standard_ATLAS_plot(input_histograms,**kwargs):

    if isinstance(input_histograms,Histogram_Wrapper):
        list_of_histograms = [input_histograms]
    elif isinstance(input_histograms,list):
        list_of_histograms = input_histograms
    else:
        print("Histogram input must be Histogram_Wrapper class or list of Histogram_Wrappers")

    normalise = kwargs["normalise"] if "normalise" in kwargs else True

    p = Single_Plot("ATLAS Plot",list_of_histograms=list_of_histograms,normalise=normalise)
    p.Initialise_Plot_Design("ATLAS")

    plt,ax,rax = p.make_plot("line-errorbar")
    p.Add_ATLAS_Label("Internal Simulation",meta_data = {"com":13,"lumi":139})

    return p,plt