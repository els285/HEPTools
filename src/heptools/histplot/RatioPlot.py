import matplotlib.pyplot as plt

from heptools.histplot.PyHist_Class import Histogram_Wrapper
from heptools.histplot.Plotting import HEP_Plot


class Ratio_Plot_ROOT(HEP_Plot):

    """
    Ratio plot class
    The mplhep style must be applied before the fig,ax,rax tuple is created
    """

    def __init__(self,plot_title,**kwargs):
        super().__init__(plot_title,**kwargs) 
        self.divisor = kwargs["divisor"] if "divisor" in kwargs else None

        assert self.divisor,"Divisor histogram not set"   

        # Construct the ratio histogram objects
        if self.Normalised:
            self.Construct_Ratio_Histograms("Norm_ROOT_hist")
        else:
            self.Construct_Ratio_Histograms("UnNorm_ROOT_hist")


    ###### Constructing list of ratio histograms

    def Construct_Ratio_Histograms(self,hist_type):   

        self.list_of_ratio_histograms = []
 
        for HW in self.list_of_histograms:
        #     if self.divisor=="self":
        #         divisor_histogram = HW
        #     else:
            divisor_histogram = self.divisor
            hist = self.Compute_Ratio(getattr(HW,hist_type),getattr(divisor_histogram,hist_type))
            wrapped=Histogram_Wrapper.Create_Wrapper(hist,HW.name,colour=HW.colour,linewidth=HW.linewidth)
            self.list_of_ratio_histograms.append(wrapped)



    @staticmethod
    def Compute_Ratio(ROOT_hist_numer,ROOT_hist_denom):

        d_hist = ROOT_hist_numer.Clone()
        d_hist.Divide(ROOT_hist_numer,ROOT_hist_denom)

        return d_hist


    def add_axis_labels(self,**kwargs):

        if "x_upper" in kwargs:
            self.axes[0].set_xlabel(kwargs["x_upper"])

        if "x_lower" in kwargs:
            self.axes[1].set_xlabel(kwargs["x_lower"])

        if "y_upper" in kwargs:
            self.axes[0].set_ylabel(kwargs["y_upper"])

        if "y_lower" in kwargs:
            self.axes[1].set_ylabel(kwargs["y_lower"])


    def Identical_Plotting(self,ax,rax,plotting_function):
        for HW in self.list_of_histograms:
            PH = HW.Norm_PyWrap_Hist if self.Normalised else HW.UnNorm_PyWrap_Hist
            plotting_function(ax,PH)

        for HW in self.list_of_ratio_histograms:
            plotting_function(rax,HW)

        # Compute extrema for ratio plot sizing
        extrema = [(min(HW.Bin_Values) , max(HW.Bin_Values)) for HW in self.list_of_ratio_histograms]

        largest_error = max(HW.Bin_Errors)

        maxY = max([a[1] for a in extrema])
        minY = min([a[0] for a in extrema])
        newYmin = (minY - largest_error) if plotting_function!=self.Step_Line else minY
        newYmax = (maxY + largest_error) if plotting_function!=self.Step_Line else maxY
        newYmin = 1 - abs(newYmin-1)*1.15
        newYmax = 1 + abs(newYmax-1)*1.15   
        # symmetric = True
        # scale = 1.25
        # if symmetric:
        #     max_value = max([abs(minY),abs(maxY)])
        #     sf = abs(1 - max_value)
        #     new_max = scale*sf


        rax.set_ylim([newYmin,newYmax])


        return ax,rax


    def One_Filled_Rest_Line(self,ax,rax,filled_hist):

        """
        Pass the overall Python wrapper as the histogram to be filled
        """

        hist4line       = [HW for HW in self.list_of_histograms if HW != filled_hist]       
        hist4line_ratio = [PY for PY in self.list_of_ratio_histograms if PY.Name != filled_hist.name]
        filled_hist_ratio_list = list(set(self.list_of_ratio_histograms)  - set(hist4line_ratio))
        assert len(filled_hist_ratio_list)==1, "Multiple histograms named the same thing: cannot identify unique histogram to plot as filled"
        filled_ratio_hist = filled_hist_ratio_list[0]

        for PH in hist4line:
            HW = PH.Norm_PyWrap_Hist if self.Normalised else PH.UnNorm_PyWrap_Hist
            ax = HEP_Plot.Line_Filled_Errors(ax,HW)

        fh_HW = filled_hist.Norm_PyWrap_Hist if self.Normalised else filled_hist.UnNorm_PyWrap_Hist
        HEP_Plot.Filled_Hist(ax,fh_HW)

        for HW in hist4line_ratio:
            rax = HEP_Plot.Line_Filled_Errors(rax,HW)
        HEP_Plot.Filled_Hist(rax,filled_ratio_hist)

        return ax,rax



    def Make_Ratio_Plot(self,plot_type,**kwargs):

        """ Assigns the correct plotting function
        Initialises the axes """

        plotting_function = self.select_plot_type(plot_type)

        """ Initialises the subplot axes and returns the fig,ax,rax tuple"""

        fig, (ax, rax) = plt.subplots(2, 1, figsize=(6,6), gridspec_kw=dict(height_ratios=[3, 1], hspace=0.1), sharex=True)
        self.fig = fig
        self.axes = (ax,rax)

        identical = False if plotting_function==self.One_Filled_Rest_Line else True

 
        if identical:
            ax,rax = self.Identical_Plotting(ax,rax,plotting_function)

        # elif not identical and plotting_function==self.One_Filled_Rest_Line and "filled" in kwargs:
        #     ax,rax = self.One_Filled_Rest_Line(ax,rax,kwargs["filled"])

        # Do the legend here
        self.do_legend()


        return plt,ax,rax





def standard_ATLAS_ratio_plot(list_of_histograms,**kwargs):

    divisor   = kwargs["divisor"]   if "divisor"   in kwargs else list_of_histograms[0]
    normalise = kwargs["normalise"] if "normalise" in kwargs else True


    p = Ratio_Plot_ROOT("A Plot",list_of_histograms=list_of_histograms,divisor=divisor,normalise=normalise)
    p.Initialise_Plot_Design("ATLAS")

    plt,ax,rax = p.Make_Ratio_Plot("line-errorbar")
    p.Add_ATLAS_Label("Internal")

    return p,plt




