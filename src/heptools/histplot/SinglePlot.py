
import matplotlib.pyplot as plt
from heptools.histplot.Plotting import HEP_Plot

class Single_Plot(HEP_Plot):
    def __init__(self,plot_title,**kwargs):
        super().__init__(plot_title,**kwargs) 


    def make_plot(self,plot_type):
        plotting_function = self.select_plot_type(plot_type)
        fig, ax = plt.subplots()
        self.axes = [ax]

        for PH in self.list_of_histograms:
            HW = PH.Norm_PyWrap_Hist if self.Normalised else PH.UnNorm_PyWrap_Hist
            plotting_function(ax,HW)
            self.do_legend()
        
        return plt, fig, ax


def standard_ATLAS_plot(list_of_histograms,**kwargs):

    normalise = kwargs["normalise"] if "normalise" in kwargs else True

    p = Single_Plot("ATLAS Plot",list_of_histograms=list_of_histograms,normalise=normalise)
    p.Initialise_Plot_Design("ATLAS")

    plt,ax,rax = p.make_plot("line-errorbar")
    p.Add_ATLAS_Label("Internal")

    return p,plt