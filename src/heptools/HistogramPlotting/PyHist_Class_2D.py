import matplotlib.pyplot as plt
import numpy as np
from ROOT import TFile,TAxis,TH1,gROOT
import os
import numpy as np
import pickle


class PyHist2D:
    """ Basic wrapper for ROOT histogram 
    Should contain no ROOT functionality, just a container for the information
    """

    def __init__(self,Name,Bin_Values,Bin_Errors,Bin_Edges,Bin_Centres,**kwargs):

        """
        This can stay in the same form as before, but meta-data probably redundant as no longer line plot
        """

        # Data
        self.Name        = Name
        self.Bin_Values  = Bin_Values  
        self.Bin_Errors  = Bin_Errors  
        self.Bin_Edges   = Bin_Edges   
        self.Bin_Centres = Bin_Centres 

        # Plot meta-data
        # self.colour       = kwargs["colour"]       if "colour"       in kwargs else "blue"
        # self.linewidth    = kwargs["linewdith"]    if "linewdith"    in kwargs else 2
        # self.legend_entry = kwargs["legend_entry"] if "legend_entry" in kwargs else ""



class Histogram_Wrapper:
    """
    Larger wrapper which contains PyHists, including methods for ROOT->numpy conversion
    Can also contain the plotting since we do not need to combine these
    """

    @staticmethod
    def get_bin_values(hist):
        outer =[]
        for binnX in range(0,hist.GetXaxis().GetNbins()):
            inner = []
            for binnY in range(0,hist.GetYaxis().GetNbins()):
                inner.append(hist.GetBinContent(binnX+1,binnY+1))
            outer.append(inner)
        arr1 = np.asarray(outer)
        arr3 = np.transpose(arr1,axes=None)
        return arr3


    @staticmethod
    def get_bin_errors(hist):
        outer =[]
        for binnX in range(0,hist.GetXaxis().GetNbins()):
            inner = []
            for binnY in range(0,hist.GetYaxis().GetNbins()):
                inner.append(hist.GetBinError(binnX+1,binnY+1))
            outer.append(inner)
        arr1 = np.asarray(outer)
        return np.flip(arr1,[0,1])


    @staticmethod
    def get_bin_edges(hist):
        x_axis = np.asarray([hist.GetXaxis().GetBinUpEdge(binn) for binn in range(0,hist.GetXaxis().GetNbins()+1)])
        y_axis = np.asarray([hist.GetYaxis().GetBinUpEdge(binn) for binn in range(0,hist.GetYaxis().GetNbins()+1)])
        return (x_axis,y_axis)


    @staticmethod
    def get_bin_centres(hist):
        x_axis = np.asarray([hist.GetXaxis().GetBinCenter(binn+1) for binn in range(0,hist.GetXaxis().GetNbins())])
        y_axis = np.asarray([hist.GetYaxis().GetBinCenter(binn+1) for binn in range(0,hist.GetYaxis().GetNbins())])
        return (x_axis,y_axis)


    @staticmethod
    def Compute_Normalised(ROOT_hist):

        """ Takes a ROOT histogram and normalised it"""

        h1dN = ROOT_hist.Clone(ROOT_hist.GetName()+"_norm")
        h1dN.Scale(1/ROOT_hist.Integral())

        return h1dN


    @staticmethod
    def normalise_by_row(ROOT_2dhist):
        sum_of_rows = []
        # Compute normalisation by row of the ROOT_2dhist
        for row in range(1, ROOT_2dhist.GetNbinsY()+1):
            sum_of_this_row = 0.
            for column in range(1, ROOT_2dhist.GetNbinsX()+1):
                sum_of_this_row += ROOT_2dhist.GetBinContent(column, row)
            sum_of_rows.append(sum_of_this_row)

        # Get normalised  output
        for row in range(1, ROOT_2dhist.GetNbinsY()+1):
            for column in range(1, ROOT_2dhist.GetNbinsX()+1):
                bin_content = 0.
                if not sum_of_rows[row-1] == 0:
                    bin_content = round((ROOT_2dhist.GetBinContent(column, row)/sum_of_rows[row-1]), 4)# * acceptance[column-1]
                ROOT_2dhist.SetBinContent(column, row, bin_content)

        return ROOT_2dhist

    def to_pandas(self):

        import pandas as pd 

        XL , YL = self.pyhist.Bin_Edges[0] , self.pyhist.Bin_Edges[1]
        Xbin_indices = [(XL[i],XL[i+1]) for i in range(0,len(XL)-1)]
        Ybin_indices = [(YL[i],YL[i+1]) for i in range(0,len(YL)-1)]

        df_og   = pd.DataFrame(self.pyhist.Bin_Values     , columns=Xbin_indices, index=Ybin_indices)
        df_norm = pd.DataFrame(self.pyhist_norm.Bin_Values, columns=Xbin_indices, index=Ybin_indices)
        return df_og,df_norm


    def plot_2d(self,normed):

        """
        Makes migration matrix ATLAS plot
        """

        import mplhep as hep
        import matplotlib.pyplot as plt

        plt.style.use(hep.style.ATLAS)

        if normed:
            pyhist = self.pyhist_norm
        else:
            pyhist = self.pyhist

        fig, ax = plt.subplots()
        X, Y = np.meshgrid(pyhist.Bin_Edges[0],pyhist.Bin_Edges[1])
        pc = ax.pcolormesh(pyhist.Bin_Edges[0],pyhist.Bin_Edges[1], pyhist.Bin_Values)#,cmap="")

        hep.atlas.text("Internal",ax=ax,loc=0)

        for i,x in enumerate(pyhist.Bin_Centres[0]):
            for j,y in enumerate(pyhist.Bin_Centres[1]):
                ax.text(y,x,pyhist.Bin_Values[i,j],ha='center', va='center',fontsize=24)

        return fig, ax, pc
    

    @classmethod
    def Create_Wrapper(self,hist,name,**kwargs):

        Bin_Values  = self.get_bin_values(hist)
        Bin_Errors  = self.get_bin_errors(hist)
        Bin_Centres = self.get_bin_centres(hist)
        Bin_Edges   = self.get_bin_edges(hist)
        return PyHist2D(name,Bin_Values,Bin_Errors,Bin_Edges,Bin_Centres,**kwargs)


    def __init__(self,ROOT_hist,name,**kwargs):

        '''
        The Histogram_Wrapper object is for converting ROOT histograms into numpy objects,
            but the data members of the class should contain no ROOT objects.
        '''

        self.ROOT_hist = ROOT_hist
        self.name      = name
        self.pyhist = self.Create_Wrapper(ROOT_hist,name)

        do_normalise_by_row = kwargs["normalise"] if "normalise" in kwargs else True

        if do_normalise_by_row:
            self.Norm_ROOT_hist = self.normalise_by_row(ROOT_hist)

        self.pyhist_norm = self.Create_Wrapper(self.Norm_ROOT_hist,name+"_normalised")







