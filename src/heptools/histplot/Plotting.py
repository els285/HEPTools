import matplotlib.pyplot as plt
import mplhep as hep
import hist
import numpy as np
from hist.intervals import ratio_uncertainty
import boost_histogram as bh
import uproot
# from matplotlib.lines import Line2D


from heptools.histplot.PyHist_Class import PyHist,Histogram_Wrapper

# from ROOT import TH1
# TH1.SetDefaultSumw2()


class HEP_Plot:

    """
    The fundamental class of HEP plot
    Includes the design specifics and the functions for each type of histogram plot
    """

    allowed_experiment_styles = ["ATLAS","ALICE","LHCb2","CMS","ATLASTex"]

    def Add_ATLAS_Label(self,label_text,**kwargs):

        assert self.axes,"Axes not defined yet. Add label after plot is made"
        ax = self.axes[0]

        if self.plot_design in ["ATLAS","ATLASTex"]:

            loc = kwargs["loc"] if "loc" in kwargs else 1

            if loc == 1:
                yaxis_limits = ax.get_ylim()
                ax.set_ylim(yaxis_limits[0],yaxis_limits[1]*1.2)
            l1 = hep.atlas.text(label_text,ax=ax,loc=loc)

            if "specific_location" in kwargs:
                #Defined from first word
                x_diff,y_diff = l1[1]._x - l1[0]._x , l1[1]._y - l1[0]._y

                l1[0]._x = kwargs["specific_location"][0]
                l1[0]._y = kwargs["specific_location"][1]
                l1[1]._x = kwargs["specific_location"][0] + x_diff
                l1[1]._y = kwargs["specific_location"][1] + y_diff

 
    def Initialise_LHC_Plot(self,experiment):

        """
        Sets the plot style based on the HEPstyle input
        """

        assert any([experiment == x for x in self.allowed_experiment_styles]), "Experiment style not defined"

        # hep.style.use(experiment)#, {'xtick.direction': 'out'}])
        if   experiment=="ATLAS": plt.style.use(hep.style.ATLAS) # This is the correct syntax for Pytohn3.6.9 version (mplhep 0.2.8)
        elif experiment=="CMS"  : plt.style.use(hep.style.CMS)
        elif experiment=="ALICE": plt.style.use(hep.style.ALICE)
        elif experiment=="LHCb2": plt.style.use(hep.style.LHCb2)
        elif experiment=="ATLASTex": 
            print("This is not supported because do not have a working TexLive distribution")
            # plt.style.use(hep.style.ATLASTex) # T


    def Initialise_Seaborn_Plot(self,**kwargs):

        import seaborn as sns
        if "style" in kwargs:
            sns.set_style(kwargs["style"])
        else:
            sns.set_style("dark")



    def Initialise_Plot_Design(self,design,**kwargs):

        self.plot_design = design

        allowed_experiment_styles = ["ATLAS","ALICE","LHCb2","CMS","ATLASTex"]

        if any([design == x for x in allowed_experiment_styles]):
            self.Initialise_LHC_Plot(design)
        elif design=="Seaborn":
            self.Initialise_Seaborn_Plot(**kwargs)



    def Add_Histograms(self,histograms2add: list):
        self.list_of_histograms = self.list_of_histograms + histograms2add




    ###################################################
    """
    Plotting functions
    """
    ###################################################


    @staticmethod
    def Step_Line(ax,PH):

        """
        Basic line histogram (which emulated mplhep.histplot)
        """

        x_binning = PH.Bin_Edges
        values    = np.concatenate((PH.Bin_Values,np.asarray([0])), axis=0) 



        ax.plot(x_binning, values,drawstyle="steps-post",color=PH.colour,label=PH.Name,linewidth=PH.linewidth)#Hist_Wrapper.linewidth)
        ax.vlines(x_binning[0],0,values[0],color=PH.colour,linewidth=PH.linewidth)#Hist_Wrapper.linewidth)

        return ax


    @staticmethod
    def Filled_Hist(ax,PH):

        """
        Basic filled histogram
        """

        x_binning = PH.Bin_Edges
        values    = np.concatenate((PH.Bin_Values,np.asarray([0])), axis=0) 

        # Required for the legend handles

        ax.plot(x_binning, values,drawstyle="steps-post",color=PH.colour,label=PH.Name,linewidth=PH.linewidth,alpha=0)#Hist_Wrapper.linewidth)
        ax.vlines(x_binning[0],0,values[0],color=PH.colour,linewidth=PH.linewidth,alpha=0)#Hist_Wrapper.linewidth)

        # The actual histogram filling

        ax.fill_between(x_binning,values,step="post", alpha=0.4,color=PH.colour)

        return ax        


    @staticmethod
    def Step_Line_Errorbar(ax,PH):

        """
        For generating a line histogram with uncapped errorbars
        """
        ax = HEP_Plot.Step_Line(ax,PH)

        ax.errorbar(PH.Bin_Centres,PH.Bin_Values,PH.Bin_Errors,elinewidth=PH.linewidth,ecolor=PH.colour,fmt='',xerr=None,linestyle='')

        return ax


    @staticmethod
    def Line_Filled_Errors(ax,PH):

        """
        For generating a line histogram with filled block error 
        """

        ax = HEP_Plot.Step_Line(ax,PH)

        for i in range(0,len(PH.Bin_Values)):

            ax.fill_between(x=[PH.Bin_Edges[i],PH.Bin_Edges[i+1]],
                            y1=PH.Bin_Values[i]+PH.Bin_Errors[i],
                            y2=PH.Bin_Values[i]-PH.Bin_Errors[i],
                            color=PH.colour,alpha=0.2)

        return ax


    def __init__(self,plot_title,**kwargs):
        self.plot_title = plot_title
        self.list_of_histograms = []

        norms = ["Normalised","normalised","Normalise","normalise"]
        self.Normalised = False
        for x in norms:
            if x in kwargs:
                self.Normalised = kwargs[x]

        self.fig  = None
        self.axes = None

        if "plot_design" in kwargs:
            self.plot_design = kwargs["plot_design"]
            self.Initalise_Plot_Design(self.plot_design,**kwargs) # Is there a way to sift a subset of kwargs i.e. **kwargs["seaborn_design_parameters"]
        if "list_of_histograms" in kwargs:
            self.list_of_histograms = kwargs["list_of_histograms"]
        if "plot_normalised" in kwargs and not kwargs["plot_normalised"]:
            self.Normalised = False


    def select_plot_type(self,plot_type):

        """Selects the type of plot from the dictionary which effectives works as a switch"""

        plot_dic = {"basic-line"        : self.Step_Line,
                    "line-errorbar"     : self.Step_Line_Errorbar,
                    "line-filled-error" : self.Line_Filled_Errors,
                    "filled-hist"       : self.Filled_Hist,
                    "one-filled-rest-line": self.One_Filled_Rest_Line
                    }

        return plot_dic[plot_type]







class Single_Plot(HEP_Plot):
    pass



class Ratio_Plot_ROOT(HEP_Plot):

    """
    Ratio plot class
    The mplhep style must be applied before the fig,ax,rax tuple is created
    """

    def __init__(self,plot_title,**kwargs):
        super().__init__(plot_title,**kwargs) 
        self.divisor = kwargs["divisor"] if "divisor" in kwargs else None
        self.Do_Ratio()



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


    ###### Constructing list of ratio histograms

    def Construct_Ratio_Histograms(self,hist_type):   

        self.list_of_ratio_histograms = []

        for HW in self.list_of_histograms:
            hist = self.Compute_Ratio(getattr(HW,hist_type),getattr(self.divisor,hist_type))
            wrapped=Histogram_Wrapper.Create_Wrapper(hist,HW.name,colour=HW.colour,linewidth=HW.linewidth)
            self.list_of_ratio_histograms.append(wrapped)



    def Do_Ratio(self):

        assert self.divisor,"Divisor histogram not set"   

        if self.Normalised:
            self.Construct_Ratio_Histograms("Norm_ROOT_hist")
        else:
            self.Construct_Ratio_Histograms("UnNorm_ROOT_hist")

    ###########################################################


    def initialise_ratio_axes(self):

        """ Initialises the subplot axes and returns the fig,ax,rax tuple"""

        fig, (ax, rax) = plt.subplots(2, 1, figsize=(6,6), gridspec_kw=dict(height_ratios=[3, 1], hspace=0.1), sharex=True)
        self.fig = fig
        self.axes = (ax,rax)
        return fig, ax, rax




    def Identical_Plotting(self,ax,rax,plotting_function):
        for PH in self.list_of_histograms:
            HW = PH.Norm_PyWrap_Hist if self.Normalised else PH.UnNorm_PyWrap_Hist
            plotting_function(ax,HW)

        for HW in self.list_of_ratio_histograms:
            plotting_function(rax,HW)

        # Compute extrema for ratio plot sizing
        extrema = [(min(HW.Bin_Values) , max(HW.Bin_Values)) for HW in self.list_of_ratio_histograms]
        # print(extrema)
        # input()
        largest_error = max(HW.Bin_Errors)
        # print(largest_error)
        # input()
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

        fig,ax,rax =self.initialise_ratio_axes()

        identical = False if plotting_function==self.One_Filled_Rest_Line else True

        # print(identical)
        # input()

        if identical:
            ax,rax = self.Identical_Plotting(ax,rax,plotting_function)

        # elif not identical and plotting_function==self.One_Filled_Rest_Line and "filled" in kwargs:
        #     ax,rax = self.One_Filled_Rest_Line(ax,rax,kwargs["filled"])

        # Do the legend here
        self.do_legend()


        return plt,ax,rax

    def do_legend(self,**kwargs):

        """
        Can redfine the legend parameters 
        """
        ax = self.axes[0]
        handles, labels = ax.get_legend_handles_labels()

        loc = kwargs["loc"] if "loc" in kwargs else "upper right"
        labels = kwargs["labels"] if "labels" in kwargs else [HW.legend_entry for HW in self.list_of_histograms]
        

        ax.legend(handles, labels, loc=loc,prop={'size': 18})



def standard_ATLAS_ratio_plot(list_of_histograms,**kwargs):

    divisor   = kwargs["divisor"]   if "divisor"   in kwargs else list_of_histograms[0]
    normalise = kwargs["normalise"] if "normalise" in kwargs else True


    p = Ratio_Plot_ROOT("A Plot",list_of_histograms=list_of_histograms,divisor=divisor,normalise=normalise)
    p.Initialise_Plot_Design("ATLAS")

    plt,ax,rax = p.Make_Ratio_Plot("line-errorbar")
    p.Add_ATLAS_Label("Internal")

    return p,plt










