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


    def __init__(self,plot_title,**kwargs):
        self.plot_title = plot_title
        self.list_of_histograms = []

        # Check for normalisation plotting setting
        norms = ["Normalised","normalised","Normalise","normalise","norm"]
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

        if any([design == x for x in self.allowed_experiment_styles]):
            self.Initialise_LHC_Plot(design)
        elif design=="Seaborn":
            self.Initialise_Seaborn_Plot(**kwargs)



    def Add_Histograms(self,histograms2add: list):
        self.list_of_histograms = self.list_of_histograms + histograms2add




    ###################################################
    ###
    ### --- Plotting functions
    ###
    ###################################################


    @staticmethod
    def data_point(ax,PH):
        x_points = PH.Bin_Centres
        y_points = PH.Bin_Values
        ax.scatter(x_points, y_points, label = 'blah',color='k')        
        ax.errorbar(x_points,y_points,yerr=PH.Bin_Errors,ls='none',color='k')

        return ax


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
    def Filled_Hist(ax,PH,**kwargs):

        """
        Basic filled histogram
        """

        x_binning = PH.Bin_Edges
        values    = np.concatenate((PH.Bin_Values,np.asarray([0])), axis=0) 
        ax.vlines(x_binning[0],0,values[0],color=PH.colour,linewidth=PH.linewidth,alpha=0)#Hist_Wrapper.linewidth)

        # Required for the legend handles

        ax.plot(x_binning, values,drawstyle="steps-post",color=PH.colour,label=PH.Name,linewidth=PH.linewidth,alpha=0)#Hist_Wrapper.linewidth)

        # The actual histogram filling
        if "y2" in kwargs:
            y2 = kwargs["y2"]
            ax.fill_between(x_binning,values,y2,step="post", alpha=0.4,color=PH.colour)
        else:
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


    def select_plot_type(self,plot_type):

        """Selects the type of plot from the dictionary which effectives works as a switch"""

        plot_dic = {"scatter-points"    : self.data_point,
                    "basic-line"        : self.Step_Line,
                    "line-errorbar"     : self.Step_Line_Errorbar,
                    "line-filled-error" : self.Line_Filled_Errors,
                    "filled-hist"       : self.Filled_Hist,
                    # "one-filled-rest-line": self.One_Filled_Rest_Line
                    }

        return plot_dic[plot_type]


    def do_legend(self,**kwargs):

        """
        Can redfine the legend parameters 
        """
        ax = self.axes[0] 
        handles, labels = ax.get_legend_handles_labels()

        loc = kwargs["loc"] if "loc" in kwargs else "upper right"
        if "labels" in kwargs:
            labels = kwargs["labels"]

        if all([hist.legend_entry is not None for hist in self.list_of_histograms]):
            labels = [hist.legend_entry for hist in self.list_of_histograms]     
     
        ax.legend(handles, labels, loc=loc,prop={'size': 18})














