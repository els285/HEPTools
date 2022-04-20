
import matplotlib.pyplot as plt
import numpy as np
from heptools.histplot.Plotting import HEP_Plot
from heptools.histplot.PyHist_Class import Histogram_Wrapper as Hist

class Data_MC_Plot(HEP_Plot):

    def __init__(self,plot_title,MC_histograms,data_histogram,**kwargs):
        super().__init__(plot_title,**kwargs) 
        self.MC_histograms = MC_histograms
        self.data_histogram = data_histogram
        self.list_of_histograms = self.MC_histograms + [self.data_histogram]

        # Check the binning
        hist0 = self.MC_histograms[0].UnNorm_PyWrap_Hist
        
        # Check for consistent bin edges
        for hist in self.MC_histograms:
            pyhist = hist.UnNorm_PyWrap_Hist
            assert np.array_equal(hist0.Bin_Edges,pyhist.Bin_Edges) , "Data-MC plot cannot be construced because MC histograms do not have coincident bin edges"

        # Stack the histograms
        list_of_arrays = [hist.UnNorm_PyWrap_Hist.Bin_Values for hist in self.MC_histograms]
        self.stacked_array = np.sum(list_of_arrays, axis=0) 


    def Data_MC_filled_hist(self,ax,PH,y1,y2):

        """
        Basic filled histogram
        """
        values = y1
        x_binning = PH.Bin_Edges

        ax.vlines(x_binning[0],0,values[0],color=PH.colour,linewidth=PH.linewidth,alpha=0)#Hist_Wrapper.linewidth)

        # Required for the legend handles
        ax.plot(x_binning, values,drawstyle="steps-post",color=PH.colour,label=PH.Name,linewidth=0,alpha=0)#Hist_Wrapper.linewidth)
        ax.fill_between(x_binning,values,y2,step="post", alpha=0.4,color=PH.colour)

        return ax    


    def make_plot(self):

        fig,ax = plt.subplots()
        if self.Normalised:
            hist_str = "Norm_PyWrap_Hist"
        else:
            hist_str = "UnNorm_PyWrap_Hist"
       
        # Plot the MC as stacked
        for i,hist in enumerate(self.MC_histograms):
            # Initialise on first iteraetion
            if i==0:
                arr = getattr(hist,hist_str).Bin_Values
                y2     = np.linspace(0,0,len(arr)+1)
                values = np.concatenate((arr,np.asarray([0])), axis=0) 

            # Drop the plot
            PH_object = getattr(hist,hist_str)
            ax = self.Data_MC_filled_hist(ax,PH_object,y1=values,y2=y2)

            # Update the paramters            
            y2 = values
            new_values    = np.concatenate((PH_object.Bin_Values,np.asarray([0])), axis=0)
            values = np.sum([y2,new_values],axis=0)
            y2 = np.concatenate((PH_object.Bin_Values,np.asarray([0])), axis=0)


        # Plot the data
        ax = self.data_point(ax,getattr(self.data_histogram,hist_str))

        self.axes = [ax]

        self.do_legend()

        return plt



def standard_ATLAS_dataMC_plot(list_of_MC_histograms,data_histogram,**kwargs):

    normalise = kwargs["normalise"] if "normalise" in kwargs else True

    p = Data_MC_Plot("ATLAS Data-MC Plot",MC_histograms=list_of_MC_histograms,data_histogram=data_histogram,normalise=normalise)
    p.Initialise_Plot_Design("ATLAS")

    plt = p.make_plot()
    p.Add_ATLAS_Label("Internal")

    return p,plt

